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
  const int ITER_TIME = 80;   // 迭代次数
  const int TRAIN_NUM = 700;  // 样本个数
  const float ALPHA = 0.015;  // 学习率
  const int NTHREAD = 4;      // 线程个数

 public:
  Logistics(const std::string &train, const std::string &predict,
            const std::string &result, const std::string &answer)
      : m_trainFile(train),
        m_predictFile(predict),
        m_resultFile(result),
        m_answerFile(answer) {}

  void LoadTrain();
  void Train();
  void Predict();
  void WriteAnswer();
  void Score();

 private:
  float Dot(const Matrix::Mat1D &mat);
  inline float sigmod(const float &z);
  void gd();
  void sgd();
  void initWeight();

 private:
  Matrix::Mat2D m_TrainData;
  std::vector<int> m_Label;
  Matrix::Mat1D m_Weight;
  std::vector<int> m_Answer;
  std::string m_trainFile;
  std::string m_predictFile;
  std::string m_resultFile;
  std::string m_answerFile;
  int m_samples;
  int m_features = 1000;
};
/***********************************************************
************************************************************
************************************************************
***************************文件处理**************************
************************************************************
************************************************************
************************************************************/

void Logistics::LoadTrain() {
  struct stat sb;
  int fd = open(m_trainFile.c_str(), O_RDONLY);
  fstat(fd, &sb);
  char *buffer = (char *)mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  Matrix::Mat1D features;
  unsigned int i = 0;
  int linesize = 6002;
  bool sign = (TRAIN_NUM != -1);
  if (sign) {
    m_TrainData = Matrix::Mat2D(TRAIN_NUM, Matrix::Mat1D(m_features));
    m_Label.resize(TRAIN_NUM);
  } else {
    int sz = sb.st_size / linesize;
    m_TrainData = Matrix::Mat2D(sz, Matrix::Mat1D(m_features));
    m_Label.resize(sz);
  }
  int x1, x2, x3, x4;
  float num;
  int pidx = 0;
  while (i < sb.st_size) {
    if (sign && pidx >= TRAIN_NUM) break;
    if (buffer[i + linesize - 1] != '\n') {
      i += linesize - 1;
      while (buffer[i] != '\n') ++i;
      ++i;
    } else {
      int idx = 0;
      while (idx < m_features) {
        x1 = buffer[i] - '0';
        x2 = buffer[i + 2] - '0';
        x3 = buffer[i + 3] - '0';
        x4 = buffer[i + 4] - '0';
        num = x1 + (float)(x2 * 100 + x3 * 10 + x4) / 1000;
        m_TrainData[pidx][idx++] = num;
        i += 6;
      }
      m_Label[pidx] = buffer[i] - '0';
      ++pidx;
      i += 2;
    }
  }
  m_samples = m_TrainData.size();
  m_features = m_TrainData[0].size();
  std::cerr << "* TrainData: (" << m_samples;
  std::cerr << ", " << m_features << ")\n";
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

  // 子线程

  auto foo = [&](int pid, int startline, int endline) {
    long long move = startline * linesize;
    Matrix::Mat1D features(m_features);
    for (int i = startline; i < endline; ++i, move += linesize) {
      int x1, x2, x3, x4;
      float num;
      int idx = 0;
      for (int j = 0; j < linesize; j += 6) {
        x1 = buffer[move + j] - '0';
        x2 = buffer[move + j + 2] - '0';
        x3 = buffer[move + j + 3] - '0';
        x4 = buffer[move + j + 4] - '0';
        num = x1 + (float)(x2 * 100 + x3 * 10 + x4) / 1000;
        features[idx++] = num;
      }
      float sigValue = this->sigmod(this->Dot(features));
      int label = (sigValue >= 0.5 ? 1 : 0);
      m_Answer[i] = label;
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

inline float Logistics::sigmod(const float &z) {
  return 1.0 / (1.0 + std::exp(-z));
}
/*
float Logistics::Dot(const Matrix::Mat1D &mat) {
  float ans = 0.0;
  for (int i = 0; i < m_features; i += 8) {
    ans += mat[i] * m_Weight[i] + mat[i + 1] * m_Weight[i + 1];
    ans += mat[i + 2] * m_Weight[i + 2] + mat[i + 3] * m_Weight[i + 3];
    ans += mat[i + 4] * m_Weight[i + 4] + mat[i + 5] * m_Weight[i + 5];
    ans += mat[i + 6] * m_Weight[i + 6] + mat[i + 7] * m_Weight[i + 7];
  }
  return ans;
}
*/
float Logistics::Dot(const Matrix::Mat1D &mat) {
  const float *p_vec1 = &mat[0];
  const float *p_vec2 = &m_Weight[0];
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
void Logistics::sgd() {
  // auto foo = [&](int cnt) {
  //   for (int i = 0; i < cnt; ++i) {
  //     int idx = Tools::RandomInt(0, m_samples - 1);
  //     auto &train = m_TrainData[idx];
  //     float sgd = this->Dot(train);
  //     float err = this->sigmod(sgd) - m_Label[idx];
  //     err *= ALPHA;
  //     for (int i = 0; i < m_features; ++i) {
  //       m_Weight[i] -= err * train[i];
  //     }
  //   }
  // };
  // std::vector<std::thread> Threads;
  // int start = 0, block = ITER_TIME / NTHREAD;
  // for (int i = 0; i < NTHREAD; ++i) {
  //   int cnt = std::min(ITER_TIME - start, block);
  //   std::thread th(foo, cnt);
  //   Threads.emplace_back(std::move(th));
  //   start += block;
  // }
  // for (auto &it : Threads) it.join();
  for (int epoch = 0; epoch < 6; ++epoch) {
    for (int i = 0; i < m_samples; ++i) {
      double rate = 5.0 / (epoch + i + 0.001) + 0.01;
      float sgd = Dot(m_TrainData[i]);
      float err = sigmod(sgd) - m_Label[i];
      for (int j = 0; j < m_features; ++j) {
        m_Weight[j] -= err * m_TrainData[i][j];
      }
    }
  }
}
void Logistics::gd() {
  int start = 0, block = m_samples / NTHREAD;
  std::vector<std::pair<int, int>> thdParam(NTHREAD);
  for (int i = 0; i < NTHREAD; ++i) {
    int end = std::min(start + block, m_samples);
    thdParam[i] = std::make_pair(start, end);
    start += block;
  }
  std::vector<std::thread> Thread(NTHREAD);
  for (int epoch = 0; epoch < ITER_TIME; ++epoch) {
    Matrix::Mat2D errAry(NTHREAD, Matrix::Mat1D(m_features, 0));
    auto foo = [&](int pid, int st, int ed) {
      for (int i = st; i < ed; ++i) {
        float sgd = Dot(m_TrainData[i]);
        sgd = sigmod(sgd) - m_Label[i];
        for (int j = 0; j < m_features; j += 8) {
          for (int k = 0; k < 8; ++k) {
            errAry[pid][j + k] += sgd * m_TrainData[i][j + k];
          }
        }
      }
    };
    int idx = 0;
    for (auto &it : thdParam) {
      Thread[idx] = std::thread(foo, idx, it.first, it.second);
      ++idx;
    }
    for (auto &it : Thread) it.join();
    float rate = 5.0 / (epoch + 1) + 0.02;
    for (auto &it : errAry) {
      for (int j = 0; j < m_features; j += 8) {
        for (int k = 0; k < 8; ++k) {
          m_Weight[j + k] -= rate * it[j + k];
        }
      }
    }
  }
}
void Logistics::initWeight() { m_Weight = Matrix::Mat1D(m_features, 1.0); }
void Logistics::Train() {
  ScopeTime t;
  initWeight();
  gd();
  // sgd();
  std::cerr << "@ train: ";
  t.LogTime();
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
  lr.LoadTrain();
  lr.Train();
  lr.Predict();
  lr.WriteAnswer();
  lr.Score();
#else
  Logistics lr(TRAIN, PREDICT, RESULT, ANSWER);
  lr.LoadTrain();
  lr.Train();
  lr.Predict();
  lr.WriteAnswer();
#endif
  std::cerr << "@ total time: ";
  t.LogTime();
  return 0;
}
