#include <fcntl.h>
#include <sys/ioctl.h>
#include <sys/ipc.h>
#include <sys/mman.h>
#include <sys/shm.h>
#include <unistd.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

struct Matrix {
  typedef std::vector<int32_t> Mat1D;
  typedef std::vector<std::vector<int32_t>> Mat2D;
};
std::ostream &operator<<(std::ostream &os, const Matrix::Mat1D &mat) {
  os << "{";
  for (size_t i = 0; i < mat.size(); ++i) {
    if (i != 0) os << ",";
    os << mat[i];
  }
  os << "}";
  return os;
}
std::ostream &operator<<(std::ostream &os, const Matrix::Mat2D &mat) {
  os << "[";
  for (auto &it : mat) {
    os << it;
  }
  os << "]";
  return os;
}

#ifdef LOCAL
const std::string TRAIN = "../data/train_data.txt";
const std::string PREDICT = "../data/test_data.txt";
const std::string RESULT = "../data/result.txt";
const std::string ANSWER = "../data/answer.txt";
#else
const std::string TRAIN = "/data/train_data.txt";
const std::string PREDICT = "/data/test_data.txt";
const std::string RESULT = "/projects/student/result.txt";
const std::string ANSWER = "/projects/student/answer.txt";
#endif

const int NTHREAD = 8;
const int SAMPLES = 5980;
const int FEATURES = 1000;
struct Node {
  int zero = 0;
  int one = 0;
  Matrix::Mat1D meansum0;
  Matrix::Mat1D meansum1;
  Node() {
    meansum0 = Matrix::Mat1D(FEATURES, 0);
    meansum1 = Matrix::Mat1D(FEATURES, 0);
  }
} ThreadData[NTHREAD];

Matrix::Mat1D MeanSum(FEATURES, 0);
Matrix::Mat1D MeanDelta(FEATURES, 0);

void HandleTrain(int pid, const char *buffer, long long start, long long end) {
  auto &data = ThreadData[pid];
  const char *ptr = buffer + start;
  long long move = 0;
  long long length = end - start;
  while (move < length) {
    const char *tail = ptr + move + 6000;
    int t = 0;
    while (*(tail + t + 1) != '\n') ++t;
    int label = *(tail + t) - '0';
    if (label) {
      ++data.one;
    } else {
      ++data.zero;
    }
    int idx = 0, num = 0;
    for (size_t i = 0; i < FEATURES; i += 8) {
      for (size_t j = 0; j < 8; ++j, move += 6, ++idx) {
        if (*(ptr + move) == '-') {
          ++move;
          num = -(*(ptr + move + 2) * 100 + *(ptr + move + 3) * 10 - 110 * '0');
        } else {
          num = *(ptr + move + 2) * 100 + *(ptr + move + 3) * 10 - 110 * '0';
        }
        if (label) {
          data.meansum1[idx] += num;
        } else {
          data.meansum0[idx] += num;
        }
      }
    }
    move += 2;
  }
}
void Train() {
  FILE *fp = fopen(RESULT.c_str(), "w");
  for (int i = 0; i < 40000; i += 2) {
    char c[2] = {' ', ' '};
    fwrite(c, 2, 1, fp);
  }
  fclose(fp);
  // 读文件
  int fd = open(TRAIN.c_str(), O_RDONLY);
  long long bufsize = SAMPLES * 6300;
  char *buffer = (char *)mmap(NULL, bufsize, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);

  long long block = bufsize / (NTHREAD + 1);
  long long start = 0, end = 0;
  std::vector<std::thread> Threads(NTHREAD);

  for (size_t i = 0; i < NTHREAD; ++i) {
    end = start + block;
    while (buffer[end] != '\n') ++end;
    ++end;
    Threads[i] = std::thread(HandleTrain, i, buffer, start, end);
    start = end;
  }
  for (auto &it : Threads) it.join();

  for (size_t i = 1; i < NTHREAD; ++i) {
    ThreadData[0].zero += ThreadData[i].zero;
    ThreadData[0].one += ThreadData[i].one;
    for (size_t j = 0; j < FEATURES; ++j) {
      ThreadData[0].meansum0[j] += ThreadData[i].meansum0[j];
      ThreadData[0].meansum1[j] += ThreadData[i].meansum1[j];
    }
  }
  int totalCount = ThreadData[0].zero + ThreadData[0].one;
  const char *ptr = buffer + end;
  int move = 0;
  for (; totalCount < SAMPLES; ++totalCount) {
    const char *tail = ptr + move + 6000;
    while (*(tail + 1) != '\n') ++tail;
    int label = *tail - '0';
    if (label) {
      ++ThreadData[0].one;
    } else {
      ++ThreadData[0].zero;
    }
    int idx = 0, num = 0;
    for (size_t i = 0; i < FEATURES; i += 8) {
      for (size_t j = 0; j < 8; ++j, move += 6, ++idx) {
        if (*(ptr + move) == '-') {
          ++move;
          num = -(*(ptr + move + 2) * 100 + *(ptr + move + 3) * 10 - 110 * '0');
        } else {
          num = *(ptr + move + 2) * 100 + *(ptr + move + 3) * 10 - 110 * '0';
        }
        if (label) {
          ThreadData[0].meansum1[idx] += num;
        } else {
          ThreadData[0].meansum0[idx] += num;
        }
      }
    }
    move += 2;
  }
  for (size_t i = 0; i < FEATURES; ++i) {
    int32_t sum0 = ThreadData[0].meansum0[i];
    int32_t sum1 = ThreadData[0].meansum1[i];
    sum0 /= ThreadData[0].zero;
    sum1 /= ThreadData[0].one;
    MeanSum[i] = sum0 + sum1;
    MeanDelta[i] = sum0 - sum1;
  }
}
void Predict(int pid) {
  long long bufsize = 20000 * 6000;
  long long block = bufsize / 8;

  int fd = open(PREDICT.c_str(), O_RDONLY);
  char *buffer = (char *)mmap(NULL, bufsize, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);

  int fd2 = open(RESULT.c_str(), O_RDWR | O_CREAT, 0666);
  char *result =
      (char *)mmap(NULL, 40000, PROT_READ | PROT_WRITE, MAP_SHARED, fd2, 0);
  close(fd2);

  buffer += (pid * block);

  int startline = pid * 2500;
  int endline = startline + 2500;

  long long move = 0;
  for (size_t line = startline; line < endline; ++line, move += 6000) {
    const char *ptr = buffer + move;
    int32_t distance = 0;
    for (size_t i = 0; i < FEATURES; ptr += 60) {
      for (size_t k = 0; k < 60; k += 6, ++i) {
        int32_t num = *(ptr + k + 2) * 200 + *(ptr + k + 3) * 20 - 220 * '0';
        int32_t sum = MeanSum[i];
        int32_t delta = MeanDelta[i];
        distance += (num - sum) * delta;
      }
    }
    char label = (distance < 0 ? '1' : '0');
    result[line << 1] = label;
    result[line << 1 | 1] = '\n';
  }
  munmap(result, 40000);
}

int main() {
  Train();

  pid_t child1 = 0, child2 = 0, child3 = 0;
  pid_t child4 = 0, child5 = 0, child6 = 0, child7 = 0;

  child1 = fork();
  if (child1) child2 = fork();
  if (child2) child3 = fork();
  if (child3) child4 = fork();
  if (child4) child5 = fork();
  if (child5) child6 = fork();
  if (child6) child7 = fork();

  int pid = 0;
  if (!child1) {
    pid = 1;
  } else if (!child2) {
    pid = 2;
  } else if (!child3) {
    pid = 3;
  } else if (!child4) {
    pid = 4;
  } else if (!child5) {
    pid = 5;
  } else if (!child6) {
    pid = 6;
  } else if (!child7) {
    pid = 7;
  } else {
    pid = 0;
  }

  if (pid <= 7) {
    Predict(pid);
    if (pid) exit(0);
  }

  // waitpid(child1, NULL, 0);
  // waitpid(child2, NULL, 0);
  // waitpid(child3, NULL, 0);
  // waitpid(child4, NULL, 0);
  // waitpid(child5, NULL, 0);
  // waitpid(child6, NULL, 0);
  // waitpid(child7, NULL, 0);

  return 0;
}
