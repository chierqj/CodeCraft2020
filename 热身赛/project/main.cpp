#include <arm_neon.h>
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

/***************************************************
 * 工具类
 * 1. RandomInt(l, r) [l, r]之间的int
 * 2. RandomDouble(l, r) [l, r]之间的double
 * 3. ShuffleVector(vector<int> vt) 随机打乱vt
 **************************************************/
class Tools {
 public:
  inline static int RandomInt(const int &l, const int &r);
  inline static float RandomDouble(const float &l, const float &r);
  inline static void ShuffleVector(std::vector<int> &vt);
};
inline int Tools::RandomInt(const int &l, const int &r) {
  static std::default_random_engine e(time(nullptr));
  std::uniform_int_distribution<int> u(l, r);
  return u(e);
}
inline float Tools::RandomDouble(const float &l, const float &r) {
  static std::default_random_engine e(time(nullptr));
  std::uniform_real_distribution<float> u(l, r);
  float ans = u(e);
  ans = floor(ans * 1000.000f + 0.5) / 1000.000f;
  return ans;
}
inline void Tools::ShuffleVector(std::vector<int> &vt) {
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::shuffle(vt.begin(), vt.end(), std::default_random_engine(seed));
}

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
 private:
  const int TRAIN_NUM = 4000;  // 样本个数
  const int NTHREAD = 4;       // 线程个数

 public:
  Logistics(const std::string &train, const std::string &predict,
            const std::string &result, const std::string &answer)
      : m_trainFile(train),
        m_predictFile(predict),
        m_resultFile(result),
        m_answerFile(answer) {}

  void Train();
  void Predict();
  void WriteAnswer();
  void Score();

 private:
  float Dot(const Matrix::Mat1D &mat1, const Matrix::Mat1D &mat2);
  void Add(Matrix::Mat1D &mat1, const Matrix::Mat1D &mat2);
  int doPredict(const Matrix::Mat1D &data);

 private:
  std::vector<int> m_Answer;
  std::string m_trainFile;
  std::string m_predictFile;
  std::string m_resultFile;
  std::string m_answerFile;
  Matrix::Mat1D m_P0Vec;
  Matrix::Mat1D m_P1Vec;
  float m_PAusuive = 0;
  float m_Log1PAusuive = 0;
  float m_Log0PAusuive = 0;
  float m_p1total = 1.0;
  float m_p0total = 1.0;
  int m_samples = 0;
  int m_positive = 0;
  int m_negative = 0;
  int m_features = 1000;
};
/***********************************************************
************************************************************
************************************************************
***************************文件处理**************************
************************************************************
************************************************************
************************************************************/

void Logistics::Train() {
  ScopeTime t;

  struct stat sb;
  int fd = open(m_trainFile.c_str(), O_RDONLY);
  fstat(fd, &sb);
  char *buffer = (char *)mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  unsigned int i = 0;
  int linesize = 6002;

  if (TRAIN_NUM != -1) {
    m_samples = TRAIN_NUM;
  } else {
    m_samples = sb.st_size / linesize;
  }

  m_P0Vec = Matrix::Mat1D(m_features, 1.0);
  m_P1Vec = Matrix::Mat1D(m_features, 1.0);

  int x1, x2, x3, x4;
  float num = 0, sum = 0;
  int pidx = 0, idx = 0;

  Matrix::Mat1D features(m_features);
  for (int i = 0; i < sb.st_size;) {
    if (pidx >= m_samples) break;
    if (buffer[i] == '-') {
      ++i;
      x1 = buffer[i] - '0';
      x2 = buffer[i + 2] - '0';
      x3 = buffer[i + 3] - '0';
      x4 = buffer[i + 4] - '0';
      num = x1 + (float)(x2 * 100 + x3 * 10 + x4) / 1000;
      features[idx++] = -num;
      sum -= num;
    } else {
      x1 = buffer[i] - '0';
      x2 = buffer[i + 2] - '0';
      x3 = buffer[i + 3] - '0';
      x4 = buffer[i + 4] - '0';
      num = x1 + (float)(x2 * 100 + x3 * 10 + x4) / 1000;
      features[idx++] = num;
      sum += num;
    }
    i += 6;
    if (idx == m_features) {
      int label = buffer[i] - '0';

      if (label == 1) {
        m_p1total += sum;
        Add(m_P1Vec, features);
        ++m_positive;
      } else {
        m_p0total += sum;
        Add(m_P0Vec, features);
        ++m_negative;
      }

      ++pidx;
      sum = 0;
      idx = 0;
      i += 2;
    }
  }

  m_PAusuive = (float)(m_positive) / (float)m_samples;
  m_Log1PAusuive = std::log(m_PAusuive);
  m_Log0PAusuive = std::log(1.0 - m_PAusuive);
  for (auto &it : m_P0Vec) it = std::log(it / m_p0total);
  for (auto &it : m_P1Vec) it = std::log(it / m_p1total);

  std::cerr << "* TrainData: (" << m_samples << ", " << m_features << ")\n";
  std::cerr << "@ train: ";
  t.LogTime();
}
void Logistics::Predict() {
  ScopeTime t;
  // 读文件
  struct stat sb;
  int fd = open(m_predictFile.c_str(), O_RDONLY);
  fstat(fd, &sb);
  char *buffer = (char *)mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);

  int linesize = 6000;
  long long linenum = sb.st_size / linesize;
  m_Answer.resize(linenum);

  float PW[10][4] = {{0.0, 0.0, 0.0, 0.0},    {1.0, 0.1, 0.01, 0.001},
                     {2.0, 0.2, 0.02, 0.002}, {3.0, 0.3, 0.03, 0.003},
                     {4.0, 0.4, 0.04, 0.004}, {5.0, 0.5, 0.05, 0.005},
                     {6.0, 0.6, 0.06, 0.006}, {7.0, 0.7, 0.07, 0.007},
                     {8.0, 0.8, 0.08, 0.008}, {9.0, 0.9, 0.09, 0.009}};
  // 子线程
  auto foo = [&](int pid, int startline, int endline) {
    long long move = startline * linesize;
    Matrix::Mat1D features(m_features);
    int x1, x2, x3, x4;
    int idx;
    float num;
    for (int i = startline; i < endline; ++i, move += linesize) {
      idx = 0;
      const char *ptr = &buffer[move];
      for (int j = 0; j < linesize; j += 6) {
        x1 = ptr[j] - '0';
        x2 = ptr[j + 2] - '0';
        x3 = ptr[j + 3] - '0';
        x4 = ptr[j + 4] - '0';
        num = PW[x1][0] + PW[x2][1] + PW[x3][2] + PW[x4][3];
        features[idx++] = num;
      }
      m_Answer[i] = this->doPredict(features);
    }
  };

  // 创建线程
  int start = 0, block = linenum / NTHREAD;
  std::vector<std::thread> Threads(NTHREAD);
  for (int i = 0; i < NTHREAD; ++i) {
    long long end = (i == NTHREAD - 1 ? linenum : start + block);
    Threads[i] = std::thread(foo, i, start, end);
    start += block;
  }
  for (auto &it : Threads) it.join();

  std::cerr << "@ predict: ";
  t.LogTime();
}

/***********************************************************
************************************************************
************************************************************
***************************训练相关**************************
************************************************************
************************************************************
************************************************************/
void Logistics::Add(Matrix::Mat1D &mat1, const Matrix::Mat1D &mat2) {
  for (int i = 0; i < m_features; ++i) {
    mat1[i] += mat2[i];
  }
}
/*
float Logistics::Dot(const Matrix::Mat1D &mat1, const Matrix::Mat1D &mat2) {
  float ans = 0.0;
  for (int i = 0; i < m_features; ++i) {
    ans += mat1[i] * mat2[i];
  }
  return ans;
}
*/
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
// void Logistics::doTrain(const Matrix::Mat1D train, const float &sum,
//                         const int label) {
//   auto add = [&](Matrix::Mat1D &mat1, const Matrix::Mat1D &mat2) {
//     for (int i = 0; i < m_features; ++i) {
//       mat1[i] += mat2[i];
//     }
//   };
//   if (label == 1) {
//     m_p1total += m_TrainSum[i];
//     add(m_P1Vec, m_TrainData[i]);
//   }

//   m_PAusuive = (float)(m_positive) / (float)m_samples;
//   m_Log1PAusuive = std::log(m_PAusuive);
//   m_Log0PAusuive = std::log(1.0 - m_PAusuive);
//   m_P0Vec = Matrix::Mat1D(m_features, 1.0);
//   m_P1Vec = Matrix::Mat1D(m_features, 1.0);
//   float p0total = 1.0, p1total = 1.0;

//   for (int i = 0; i < m_samples; ++i) {
//     if (m_Label[i] == 1) {
//     } else {
//       p0total += m_TrainSum[i];
//       add(m_P0Vec, m_TrainData[i]);
//     }
//   }

//   for (auto &it : m_P0Vec) it = std::log(it / p0total);
//   for (auto &it : m_P1Vec) it = std::log(it / p1total);
// }
int Logistics::doPredict(const Matrix::Mat1D &data) {
  // float sigValue = this->sigmod(this->Dot(features));
  // int label = (sigValue >= 0.5 ? 1 : 0);

  float p1 = Dot(data, m_P1Vec) + m_Log1PAusuive;
  float p0 = Dot(data, m_P0Vec) + m_Log0PAusuive;
  return (p1 > p0 ? 1 : 0);
}

void Logistics::WriteAnswer() {
  ScopeTime t;
  FILE *fp = fopen(m_resultFile.c_str(), "w");
  for (auto &label : m_Answer) {
    char c[2];
    c[0] = label + '0';
    c[1] = '\n';
    fwrite(c, 2, 1, fp);
  }
  fclose(fp);
  std::cerr << "@ write result: ";
  t.LogTime();
}
void Logistics::Score() {
  std::ifstream fin(m_answerFile);
  assert(fin);
  int x, index = 0;
  float ac = 0, tol = 0;
  while (fin >> x) {
    if (x == m_Answer[index++]) {
      ++ac;
    }
    ++tol;
  }
  fin.close();
  std::cerr << "@ threads: " << NTHREAD << "\n";
  std::cerr << "@ accuracy: " << ac * 100 / tol << "%\n";
}

int main() {
  std::cerr << std::fixed << std::setprecision(3);
  const std::string TRAIN = "/data/train_data.txt";
  const std::string PREDICT = "/data/test_data.txt";
  const std::string RESULT = "/projects/student/result.txt";
  const std::string ANSWER = "/projects/student/answer.txt";
  const std::string LOCAL_TRAIN = "../data/train_data.txt";
  const std::string LOCAL_PREDICT = "../data/test_data.txt";
  const std::string LOCAL_RESULT = "../data/result.txt";
  const std::string LOCAL_ANSWER = "../data/answer.txt";

  ScopeTime t;

#ifdef LOCAL
  Logistics lr(LOCAL_TRAIN, LOCAL_PREDICT, LOCAL_RESULT, LOCAL_ANSWER);
  lr.Train();
  lr.Predict();
  lr.WriteAnswer();
  lr.Score();
#else
  Logistics lr(TRAIN, PREDICT, RESULT, ANSWER);
  lr.Train();
  lr.Predict();
  lr.WriteAnswer();
#endif
  std::cerr << "@ total time: ";
  t.LogTime();
  return 0;
}
