#include <fcntl.h>
#include <sys/ioctl.h>
#include <sys/ipc.h>
#include <sys/mman.h>
#include <sys/shm.h>
#include <sys/wait.h>
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

#ifdef LOCAL
#define TRAIN "../data/train_data.txt"
#define PREDICT "../data/test_data.txt"
#define RESULT "../data/result.txt"
#define ANSWER "../data/answer.txt"
#else
#define TRAIN "/data/train_data.txt"
#define PREDICT "/data/test_data.txt"
#define RESULT "/projects/student/result.txt"
#define ANSWER "/projects/student/answer.txt"
#endif

const int NTHREAD = 8;
const int SAMPLES = 1996;  // 6966(6068)
const int FEATURES = 1000;
struct Node {
  int zero = 0;
  int one = 0;
  int32_t meansum0[FEATURES] = {0};
  int32_t meansum1[FEATURES] = {0};
} ThreadData[NTHREAD];

int32_t MeanSum[FEATURES] = {0};
int32_t MeanDelta[FEATURES] = {0};

void Init() {
  FILE *fp = fopen(RESULT, "w");
  for (int i = 0; i < 40000; i += 2) {
    char c[2] = {' ', ' '};
    fwrite(c, 2, 1, fp);
  }
  fclose(fp);
}

inline int32_t GetNumber(const char *ptr, int pos) {
  // return (*(ptr + pos + 2) - '0') * 100;
  return (*(ptr + pos + 2) - '0') * 100 + (*(ptr + pos + 3) - '0') * 10;
}
void HandleTrain(int pid, const char *buffer, long long start, long long end) {
  auto &data = ThreadData[pid];
  const char *ptr = buffer + start;
  while (start < end) {
    const char *tail = ptr + 6000;
    while (*(tail + 1) != '\n') ++tail;
    bool sign = (*tail == '1' ? true : false);
    if (sign) {
      ++data.one;
    } else {
      ++data.zero;
    }

    int32_t num = 0;
    for (size_t i = 0; i < FEATURES;) {
      int32_t move = 0;
      for (size_t j = 0, pos = 0; j < 8; ++j, ++i, pos += 6) {
        if (*(ptr + pos) == '-') {
          ++move;
          ++pos;
          num = -GetNumber(ptr, pos);
        } else {
          num = GetNumber(ptr, pos);
        }
        if (sign) {
          data.meansum1[i] += num;
        } else {
          data.meansum0[i] += num;
        }
      }
      ptr += (48 + move);
      start += (48 + move);
    }
    ptr += 2;
    start += 2;
  }
}
void Train() {
  // 读文件
  int fd = open(TRAIN, O_RDONLY);
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
  int32_t num = 0;
  for (; totalCount < SAMPLES; ++totalCount) {
    const char *tail = ptr + 6000;
    while (*(tail + 1) != '\n') ++tail;
    bool sign = (*tail == '1' ? true : false);
    if (sign) {
      ++ThreadData[0].one;
    } else {
      ++ThreadData[0].zero;
    }
    int32_t num = 0;
    for (size_t i = 0; i < FEATURES;) {
      int32_t move = 0;
      for (size_t j = 0, pos = 0; j < 8; ++j, ++i, pos += 6) {
        if (*(ptr + pos) == '-') {
          ++move;
          ++pos;
          num = -GetNumber(ptr, pos);
        } else {
          num = GetNumber(ptr, pos);
        }
        if (sign) {
          ThreadData[0].meansum1[i] += num;
        } else {
          ThreadData[0].meansum0[i] += num;
        }
      }
      ptr += (48 + move);
    }
    ptr += 2;
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
  long long bufsize = 120000000;
  long long block = 15000000;

  int fd = open(PREDICT, O_RDONLY);
  char *buffer = (char *)mmap(NULL, bufsize, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);

  int fd2 = open(RESULT, O_RDWR | O_CREAT, 0666);
  char *result =
      (char *)mmap(NULL, 40000, PROT_READ | PROT_WRITE, MAP_SHARED, fd2, 0);
  close(fd2);

  buffer += (pid * block);

  int startline = pid * 2500;
  int endline = startline + 2500;

  const char *ptr = buffer;
  int32_t distance = 0, num;

  for (size_t line = startline; line < endline; ++line) {
    if (ptr[2] >= '2') {
      result[line << 1] = '1';
      result[line << 1 | 1] = '\n';
      ptr += 6000;
      continue;
    }
    distance = 0;
    for (size_t i = 0; i < FEATURES; ptr += 60) {
      for (size_t k = 0; k < 60; k += 6, ++i) {
        // num = (*(ptr + k + 2) - '0') * 200;
        num = (*(ptr + k + 2) - '0') * 200 + (*(ptr + k + 3) - '0') * 20;
        distance += (num - MeanSum[i]) * MeanDelta[i];
      }
    }
    result[line << 1] = (distance < 0 ? '1' : '0');
    result[line << 1 | 1] = '\n';
  }
  // munmap(result, 40000);
}

int main() {
  Init();
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
    exit(0);
  }
  /*
  waitpid(child1, NULL, 0);
  waitpid(child2, NULL, 0);
  waitpid(child3, NULL, 0);
  waitpid(child4, NULL, 0);
  waitpid(child5, NULL, 0);
  waitpid(child6, NULL, 0);
  waitpid(child7, NULL, 0);
  */
  return 0;
}
