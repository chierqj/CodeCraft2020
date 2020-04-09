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
#include <unordered_map>
#include <vector>
#include <arm_neon.h>


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

  void InitData();
  void Train();
  void Predict();
  void WriteAnswer();
  void Score();
  float Dot(const Matrix::Mat1D &mat1, const Matrix::Mat1D &mat2);
  int doPredict(const Matrix::Mat1D &data);
  inline float getNumber(const char *ptr, int pos);

 public:
  int m_samples = 1580;
  int m_features = 400;
  const int NTHREAD = 8;  // 线程个数

 private:
  Matrix::Mat1D m_MeanLabel[2];
  int m_CountLabel[2];
  float m_PredictP0 = 0;
  float m_PredictP1 = 0;

  Matrix::Mat1D m_PredictSum;
  Matrix::Mat1D m_PredictDelta;
  Matrix::Mat1D m_Answer;

 private:
  std::string m_trainFile;
  std::string m_predictFile;
  std::string m_resultFile;
  std::string m_answerFile;
};

void Logistics::InitData() {
  m_MeanLabel[0] = Matrix::Mat1D(m_features, 1.0);
  m_MeanLabel[1] = Matrix::Mat1D(m_features, 1.0);
  m_CountLabel[0] = 0;
  m_CountLabel[1] = 0;

  m_PredictSum = Matrix::Mat1D(20000, 0);
  m_PredictDelta = Matrix::Mat1D(20000, 0);
  m_Answer = Matrix::Mat1D(20000, 0);
}
inline float Logistics::getNumber(const char *ptr, int pos) {
  float num = (ptr[pos] - '0') + (ptr[pos + 2] - '0') * 0.1 +
              (ptr[pos + 3] - '0') * 0.01;
  return num;
}

void Logistics::Train() {
  ScopeTime t;

  struct Node {
    Matrix::Mat1D featureSum[2];
    int labelCount[2];
    float sum[2];
    Node(int n) {
      featureSum[0] = Matrix::Mat1D(n, 0);
      featureSum[1] = Matrix::Mat1D(n, 0);
      labelCount[0] = 0;
      labelCount[1] = 0;
      sum[0] = 0;
      sum[1] = 0;
    }
  };
  float total[2] = {1.0, 1.0};

  std::vector<Node> ThreadData(NTHREAD, Node(m_features));

  auto foo = [&](int pid, long long start, long long end) {
    std::ifstream thfin(m_trainFile);
    thfin.seekg(start);
    std::string line;
    auto &thData = ThreadData[pid];
    int label = 0, pos = 0;
    while (start < end && std::getline(thfin, line)) {
      label = line.back() - '0';
      ++thData.labelCount[label];
      pos = 0;
      float tmpsum = 0;
      // for (int i = 0; i < m_features; ++i, pos += 6) {
      //   bool sign = false;
      //   if (line[pos] == '-') {
      //     ++pos;
      //     sign = true;
      //   }
      //   float num = getNumber(line, pos);
      //   if (sign) num = -num;
      //   thData.featureSum[label][i] += num;
      //   tmpsum += num;
      // }

      for (int i = 0; i < m_features; pos += 48) {
        const char *ptr = &line[pos];
        for (int k = 0, move = 0; k < 8; ++k, ++i, move += 6) {
          bool sign = false;
          if (ptr[move] == '-') {
            ++pos;
            ++move;
            sign = true;
          }
          float num = getNumber(ptr, move);
          if (sign) num = -num;
          thData.featureSum[label][i] += num;
          tmpsum += num;
        }
      }

      start += line.size() + 1;
      thData.sum[label] += tmpsum;
    }
  };

  std::ifstream fin(m_trainFile);
  long long bufsize = m_samples * 6500;
  long long block = bufsize / (NTHREAD + 1);
  long long start = 0, end = 0;
  std::vector<std::thread> Threads(NTHREAD);

  for (int i = 0; i < NTHREAD; ++i) {
    end = start + block;
    fin.seekg(end, std::ios::beg);
    std::string line;
    std::getline(fin, line);
    end = fin.tellg();
    Threads[i] = std::thread(foo, i, start, end);
    start = end;
  }
  for (auto &it : Threads) it.join();

  for (int i = 0; i < NTHREAD; ++i) {
    const auto &thData = ThreadData[i];
    m_CountLabel[0] += thData.labelCount[0];
    m_CountLabel[1] += thData.labelCount[1];
    total[0] += thData.sum[0];
    total[1] += thData.sum[1];
    for (int j = 0; j < m_features; ++j) {
      m_MeanLabel[0][j] += thData.featureSum[0][j];
      m_MeanLabel[1][j] += thData.featureSum[1][j];
    }
  }
  int totalCount = m_CountLabel[0] + m_CountLabel[1];
  for (; totalCount < m_samples; ++totalCount) {
    std::string line;
    std::getline(fin, line);
    int label = line.back() - '0';
    ++m_CountLabel[label];
    int pos = 0;
    float tmpsum = 0;
    // for (int i = 0; i < m_features; ++i, pos += 6) {
    //   bool sign = false;
    //   if (line[pos] == '-') {
    //     ++pos;
    //     sign = true;
    //   }
    //   float num = getNumber(line, pos);
    //   if (sign) num = -num;
    //   m_MeanLabel[label][i] += num;
    //   tmpsum += num;
    // }

    for (int i = 0; i < m_features; pos += 48) {
      const char *ptr = &line[pos];
      for (int k = 0, move = 0; k < 8; ++k, ++i, move += 6) {
        bool sign = false;
        if (ptr[move] == '-') {
          ++pos;
          ++move;
          sign = true;
        }
        float num = getNumber(ptr, move);
        if (sign) num = -num;
        m_MeanLabel[label][i] += num;
        tmpsum += num;
      }
    }
    total[label] += tmpsum;
  }
  fin.close();

  m_PredictP0 = (float)m_CountLabel[0] / (float)totalCount;
  m_PredictP1 = 1.0 - m_PredictP0;
  m_PredictP0 = std::log(m_PredictP0);
  m_PredictP1 = std::log(m_PredictP1);

  for (int i = 0; i < m_features; ++i) {
    m_MeanLabel[0][i] = std::log(m_MeanLabel[0][i] / total[0]);
    m_MeanLabel[1][i] = std::log(m_MeanLabel[1][i] / total[1]);
  }

  std::cerr << "@ train: ";
  t.LogTime();
}
// float Logistics::Dot(const Matrix::Mat1D &mat1, const Matrix::Mat1D &mat2) {
//   float ans = 0.0;
//   for (int i = 0; i < m_features; ++i) {
//     ans += mat1[i] * mat2[i];
//   }
//   return ans;
// }
float Logistics::Dot(const Matrix::Mat1D &mat1, const Matrix::Mat1D &mat2) {
  const float *p_vec1 = &mat1[0];
  const float *p_vec2 = &mat2[0];
  float sum = 0;
  float32x4_t sum_vec = vdupq_n_f32(0), left_vec, right_vec;
  for (int i = 0; i < m_features; i += 4) {
    left_vec = vld1q_f32(p_vec1 + i);
    right_vec = vld1q_f32(p_vec2 + i);
    sum_vec = vmlaq_f32(sum_vec, left_vec, right_vec);
  }
  float32x2_t r = vadd_f32(vget_high_f32(sum_vec), vget_low_f32(sum_vec));
  sum += vget_lane_f32(vpadd_f32(r, r), 0);
  return sum;
}

int Logistics::doPredict(const Matrix::Mat1D &data) {
  float p1 = Dot(data, m_MeanLabel[1]) + m_PredictP1;
  float p0 = Dot(data, m_MeanLabel[0]) + m_PredictP0;
  return (p1 > p0 ? 1 : 0);
}
// void Logistics::Predict() {
//   ScopeTime t;

//   auto foo = [&](int start, int end) {
//     std::ifstream thfin(m_predictFile);
//     thfin.seekg((long long)start * 6000L, std::ios::beg);
//     std::string line;
//     Matrix::Mat1D features(m_features);
//     while (start < end && std::getline(thfin, line)) {
//       int pos = 0;
//       // for (int i = 0; i < m_features; ++i, pos += 6) {
//       //   features[i] = getNumber(line, pos);
//       // }
//       for (int i = 0; i < m_features; pos += 60) {
//         const char *ptr = &line[pos];
//         for (int k = 0; k < 60; k += 6, ++i) {
//           features[i] = getNumber(ptr, k);
//           // x1 = ptr[k + 2] - '0';
//           // x2 = ptr[k + 3] - '0';
//           // int num = x1 * 200 + x2 * 20;
//           // distance += (num - m_PredictSum[i]) * m_PredictDelta[i];
//         }
//       }

//       m_Answer[start++] = doPredict(features);
//     }
//   };

//   std::ifstream fin(m_predictFile);
//   long long bufsize = fin.seekg(0, std::ios::end).tellg();
//   int linecount = bufsize / 6000;
//   int block = linecount / NTHREAD;
//   int start = 0, end = block;
//   std::vector<std::thread> Threads(NTHREAD);
//   for (int i = 0; i < NTHREAD; ++i) {
//     end = (i == NTHREAD - 1 ? linecount : start + block);
//     Threads[i] = std::thread(foo, start, end);
//     start += block;
//   }
//   for (auto &it : Threads) it.join();

//   fin.close();

//   std::cerr << "@ predict: ";
//   t.LogTime();
// }

void Logistics::Predict() {
  ScopeTime t;
  struct stat sb;
  int fd = open(m_predictFile.c_str(), O_RDONLY);
  fstat(fd, &sb);
  char *buffer = (char *)mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);

  int linesize = 6000;
  long long linenum = sb.st_size / linesize;
  auto foo = [&](int startline, int endline) {
    long long move = startline * linesize;
    Matrix::Mat1D features(m_features);
    for (int line = startline; line < endline; ++line, move += linesize) {
      int distance = 0;
      for (int i = 0, pos = 0; i < m_features; pos += 60) {
        const char *ptr = &buffer[move + pos];
        for (int k = 0; k < 60; k += 6, ++i) {
          features[i] = getNumber(ptr, k);
        }
      }
      m_Answer[line] = doPredict(features);
    }
  };
  int start = 0, block = linenum / NTHREAD;
  std::vector<std::thread> Threads(NTHREAD);
  for (int i = 0; i < NTHREAD; ++i) {
    long long end = (i == NTHREAD - 1 ? linenum : start + block);
    Threads[i] = std::thread(foo, start, end);
    start += block;
  }
  for (auto &it : Threads) it.join();

  std::cerr << "@ predict: ";
  t.LogTime();
}

void Logistics::WriteAnswer() {
  FILE *fp = fopen(m_resultFile.c_str(), "w");
  for (auto &label : m_Answer) {
    char c[2];
    c[0] = label + '0';
    c[1] = '\n';
    fwrite(c, 2, 1, fp);
  }
  fclose(fp);
}
void Logistics::Score() {
  std::ifstream fin(m_answerFile);
  int x, index = 0;
  float ac = 0, tol = 0;
  while (fin >> x) {
    if (x == m_Answer[index++]) {
      ++ac;
    }
    ++tol;
  }
  fin.close();
  std::cerr << "@ accuracy: " << ac * 100 / tol << "%\n";
}

int main() {
  std::cerr << std::fixed << std::setprecision(3);
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(nullptr);

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
  lr.InitData();
  lr.Train();
  lr.Predict();
  lr.WriteAnswer();
  lr.Score();

#else
  Logistics lr(TRAIN, PREDICT, RESULT, ANSWER);
  lr.InitData();
  lr.Train();
  lr.Predict();
  lr.WriteAnswer();
#endif
  return 0;
}
