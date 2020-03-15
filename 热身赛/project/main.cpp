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
  inline static double RandomDouble(const double &l, const double &r);
  inline static void ShuffleVector(std::vector<int> &vt);
};
inline int Tools::RandomInt(const int &l, const int &r) {
  static std::default_random_engine e(time(nullptr));
  std::uniform_int_distribution<int> u(l, r);
  return u(e);
}
inline double Tools::RandomDouble(const double &l, const double &r) {
  static std::default_random_engine e(time(nullptr));
  std::uniform_real_distribution<double> u(l, r);
  double ans = u(e);
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
  const int ITER_TIME = 50000;  // 迭代次数
  const int TRAIN_NUM = 5000;   // 样本个数
  const double ALPHA = 0.015;   // 学习率
  const int NTHREAD = 4;        // 线程个数

 public:
  Logistics(const std::string &train, const std::string &predict,
            const std::string &result, const std::string &answer)
      : m_trainFile(train),
        m_predictFile(predict),
        m_resultFile(result),
        m_answerFile(answer) {}

  inline void LoadData();
  void Train();
  void Predict();
  void Score();

 private:
  void loadTrain();
  void loadPredict();
  double Dot(const Matrix::Mat1D &mat);
  inline double sigmod(const double &z);
  void handleThread(const char *buffer, int pid, unsigned int start,
                    unsigned int end, unsigned int linesize);

 private:
  Matrix::Mat2D m_TrainData;
  std::vector<Matrix::Mat2D> m_ThreadData;
  std::vector<int> m_Label;
  Matrix::Mat2D m_PredictData;
  Matrix::Mat1D m_Weight;
  std::vector<int> m_Answer;
  std::string m_trainFile;
  std::string m_predictFile;
  std::string m_resultFile;
  std::string m_answerFile;
  int m_samples;
  int m_features;
};
void Logistics::loadTrain() {
  struct stat sb;
  int fd = open(m_trainFile.c_str(), O_RDONLY);
  fstat(fd, &sb);
  char *buffer = (char *)mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  Matrix::Mat1D features;
  unsigned int i = 0;
  int cnt = 0;
  bool sign = (TRAIN_NUM != -1);
  int x1, x2, x3, x4;
  double num;
  while (i < sb.st_size) {
    if (sign && cnt >= TRAIN_NUM) break;
    if (buffer[i] == '-') {
      features.clear();
      while (buffer[i] != '\n') ++i;
      ++i;
    } else if (buffer[i + 1] == '\n') {
      int x = buffer[i] - '0';
      m_Label.emplace_back(x);
      m_TrainData.emplace_back(features);
      features.clear();
      ++cnt;
      i += 2;
    } else {
      x1 = buffer[i] - '0';
      x2 = buffer[i + 2] - '0';
      x3 = buffer[i + 3] - '0';
      x4 = buffer[i + 4] - '0';
      num = x1 + (double)(x2 * 100 + x3 * 10 + x2) / 1000;
      features.emplace_back(num);
      i += 6;
    }
  }
  m_samples = m_TrainData.size();
  m_features = m_TrainData[0].size();
  std::cerr << "* TrainData: (" << m_samples;
  std::cerr << ", " << m_features << ")\n";
}
void Logistics::handleThread(const char *buffer, int pid, unsigned int start,
                             unsigned int end, unsigned int linesize) {
  for (unsigned i = start; i < end; ++i) {
    int x1, x2, x3, x4;
    double num;
    Matrix::Mat1D features;
    unsigned int index = i * linesize;
    for (int j = 0; j < linesize; j += 6) {
      x1 = buffer[index + j] - '0';
      x2 = buffer[index + j + 2] - '0';
      x3 = buffer[index + j + 3] - '0';
      x4 = buffer[index + j + 4] - '0';
      num = x1 + (double)(x2 * 100 + x3 * 10 + x2) / 1000;
      features.emplace_back(num);
    }
    m_ThreadData[pid].emplace_back(features);
  }
}
void Logistics::loadPredict() {
  struct stat sb;
  int fd = open(m_predictFile.c_str(), O_RDONLY);
  fstat(fd, &sb);
  char *buffer = (char *)mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  Matrix::Mat1D features;
  int p = 0, linesize = 0;
  int x1, x2, x3, x4;
  double num;
  while (p < sb.st_size) {
    x1 = buffer[p] - '0';
    x2 = buffer[p + 2] - '0';
    x3 = buffer[p + 3] - '0';
    x4 = buffer[p + 4] - '0';
    num = x1 + (double)(x2 * 100 + x3 * 10 + x2) / 1000;
    features.emplace_back(num);
    if (buffer[p + 5] == '\n') {
      m_PredictData.emplace_back(features);
      linesize = p + 6;
      features.clear();
      break;
    }
    p += 6;
  }
  int totalLine = sb.st_size / linesize;
  int block = totalLine / NTHREAD;
  int start = 1;
  std::vector<std::thread> Threads;
  m_ThreadData.resize(NTHREAD);
  for (int i = 0; i < NTHREAD; ++i) {
    int end = (i == NTHREAD - 1 ? totalLine : start + block);
    std::thread th(&Logistics::handleThread, this, buffer, i, start, end,
                   linesize);
    Threads.emplace_back(std::move(th));
    start += block;
  }
  for (auto &it : Threads) it.join();
  for (auto &it : m_ThreadData) {
    m_PredictData.insert(m_PredictData.end(), it.begin(), it.end());
  }
  std::cerr << "* PredictData: (" << m_PredictData.size();
  std::cerr << ", " << m_PredictData[0].size() << ")\n";
}
inline void Logistics::LoadData() {
  this->loadTrain();
  this->loadPredict();
}
inline double Logistics::sigmod(const double &z) {
  return 1.0 / (1.0 + std::exp(-z));
}
double Logistics::Dot(const Matrix::Mat1D &mat) {
  double ans = 0.0L;
  int i = 0;
  while (i < m_features - 4) {
    ans += mat[i] * m_Weight[i] + mat[i + 1] * m_Weight[i + 1] +
           mat[i + 2] * m_Weight[i + 2] + mat[i + 3] * m_Weight[i + 3];
    i += 4;
  }
  while (i < m_features) {
    ans += mat[i] * m_Weight[i];
    ++i;
  }
  return ans;
}
void Logistics::Train() {
  m_Weight = Matrix::Mat1D(m_features, 1.0);
  for (int epoch = 0; epoch < ITER_TIME; ++epoch) {
    int index = Tools::RandomInt(0, m_samples - 1);
    auto &train = m_TrainData[index];
    double sgd = this->Dot(train);
    double err = this->sigmod(sgd) - m_Label[index];
    err *= ALPHA;

    for (int i = 0; i < m_features; ++i) {
      m_Weight[i] -= err * train[i];
    }
  }
}

void Logistics::Predict() {
  auto getLabel = [&](const Matrix::Mat1D &data) {
    double sigValue = this->sigmod(this->Dot(data));
    return (sigValue >= 0.5 ? 1 : 0);
  };
  FILE *fp = fopen(m_resultFile.c_str(), "w");
  for (auto &test : m_PredictData) {
    int label = getLabel(test);
    m_Answer.emplace_back(label);
    char c[2];
    c[0] = label + '0';
    c[1] = '\n';
    fwrite(c, 2, 1, fp);
  }
  fclose(fp);
}
void Logistics::Score() {
  std::ifstream fin(m_answerFile);
  assert(fin);
  int x, index = 0;
  double ac = 0, tol = 0;
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
  lr.LoadData();
  lr.Train();
  lr.Predict();
  lr.Score();
#else
  Logistics lr(TRAIN, PREDICT, RESULT, ANSWER);
  lr.LoadData();
  lr.Train();
  lr.Predict();
#endif
  std::cerr << "@ total time: ";
  t.LogTime();
  return 0;
}
