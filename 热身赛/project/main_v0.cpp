#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

class ScopeTime {
 public:
  ScopeTime() : m_begin(std::chrono::high_resolution_clock::now()) {}
  void LogTime() const {
    auto t = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - m_begin);
    float elapsed = (float)(t.count() * 1.0) / 1000.0;
    std::cerr << elapsed << "s\n";
  }

 private:
  std::chrono::time_point<std::chrono::high_resolution_clock> m_begin;
};

struct Matrix {
  typedef std::vector<float> Mat1D;
  typedef std::vector<std::vector<float>> Mat2D;
};
std::ostream &operator<<(std::ostream &os, const Matrix::Mat1D &mat) {
  os << "{";
  for (int i = 0; i < mat.size(); ++i) {
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

class Logistics {
 public:
  Logistics(const std::string &train, const std::string &predict,
            const std::string &result, const std::string &answer)
      : m_trainFile(train),
        m_predictFile(predict),
        m_resultFile(result),
        m_answerFile(answer) {}

  void Run();
  void LoadData();
  // void Predict();
  // void Score();

 private:
  Matrix::Mat2D m_TrainData;
  std::string m_trainFile;
  std::string m_predictFile;
  std::string m_resultFile;
  std::string m_answerFile;
};
void Logistics::Run() { LoadData(); }
void Logistics::LoadData() {
  struct stat sb;
  int fd = open(m_predictFile.c_str(), O_RDONLY);
  fstat(fd, &sb); /* 取得文件大小 */
  char *data;
  std::cerr << sb.st_size << "\n";
  data = (char *)mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);

  bool negative = false, dot = false;
  float integer = 0, decimal = 0, multiple = 10;

  auto init = [&]() {
    negative = dot = false;
    integer = decimal = 0;
    multiple = 10;
  };
  auto number = [&]() {
    float num = integer + decimal;
    if (negative) num = -num;
    if (num < 0) num = 0;
    if (num > 1) num = 1;
    return num;
  };

  Matrix::Mat1D features;
  for (long long i = 0; i < sb.st_size; i += 6) {
    float num = data[i] - '0' + data[i + 2] / 10 + data[i + 3] / 100 +
                data[i + 4] / 1000;
    features.emplace_back(num);
    if (data[i + 5] == '\n') {
      m_TrainData.emplace_back(features);
      features.clear();
    }
  }
  std::cerr << m_TrainData.size() << ", " << m_TrainData[0].size() << "\n";
  // int cnt = 0;
  // bool negative = false, dot = false;
  // float integer = 0, decimal = 0, multiple = 10;
  // Matrix::Mat1D features;
  // auto init = [&]() {
  //   negative = dot = false;
  //   integer = decimal = 0;
  //   multiple = 10;
  // };
  // auto number = [&]() {
  //   float num = integer + decimal;
  //   if (negative) num = -num;
  //   if (num < 0) num = 0;
  //   if (num > 1) num = 1;
  //   return num;
  // };
  // for (long long i = 0; i < size; ++i) {
  //   if (cnt >= skip) break;
  //   char v = data[i];
  //   if (v == '\n') {
  //     features.emplace_back(number());
  //     m_TrainData.emplace_back(features);
  //     ++cnt;
  //     features.clear();
  //     init();
  //   } else if (v == ',') {
  //     features.emplace_back(number());
  //     init();
  //   } else if (v == '-') {
  //     negative = true;
  //   } else if (v == '.') {
  //     dot = true;
  //   } else {
  //     float x = v - '0';
  //     if (!dot) {
  //       integer = integer * 10 + x;
  //     } else {
  //       decimal = decimal + x / multiple;
  //       multiple *= 10;
  //     }
  //   }
  // }
}

int main() {
  ScopeTime t;
  const std::string TRAIN = "/data/train_data.txt";
  const std::string PREDICT = "/data/test_data.txt";
  const std::string RESULT = "/projects/student/result.txt";
  const std::string ANSWER = "/projects/student/answer.txt";
  const std::string LOCAL_TRAIN = "../data/train_data.txt";
  const std::string LOCAL_PREDICT = "../data/test_data.txt";
  const std::string LOCAL_RESULT = "../data/result.txt";
  const std::string LOCAL_ANSWER = "../data/answer.txt";
#ifdef LOCAL
  Logistics lr(LOCAL_TRAIN, LOCAL_PREDICT, LOCAL_RESULT, LOCAL_ANSWER);
#else
  Logistics lr(TRAIN, PREDICT, RESULT, ANSWER);
#endif
  lr.Run();
  t.LogTime();
  return 0;
}