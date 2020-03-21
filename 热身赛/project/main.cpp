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
#include <unordered_map>
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
 private:
  const float m_PW[10][4] = {{0.0, 0.0, 0.0, 0.0},    {1.0, 0.1, 0.01, 0.001},
                             {2.0, 0.2, 0.02, 0.002}, {3.0, 0.3, 0.03, 0.003},
                             {4.0, 0.4, 0.04, 0.004}, {5.0, 0.5, 0.05, 0.005},
                             {6.0, 0.6, 0.06, 0.006}, {7.0, 0.7, 0.07, 0.007},
                             {8.0, 0.8, 0.08, 0.008}, {9.0, 0.9, 0.09, 0.009}};

 private:
  std::vector<int> m_Answer;
  std::string m_trainFile;
  std::string m_predictFile;
  std::string m_resultFile;
  std::string m_answerFile;
  int m_samples = 0;
  int m_features = 400;
  const int TRAIN_NUM = 1580;  // 样本个数
  const int NTHREAD = 4;       // 线程个数

  struct Node {
    Matrix::Mat2D train;
    Matrix::Mat1D sum;
    std::vector<int> label;
    int count = 0;
    Node(int n, int m) {
      train = Matrix::Mat2D(n, Matrix::Mat1D(m));
      label.resize(n);
      sum.resize(n);
    }
  };

  std::vector<Node> m_ThreadData;
  Matrix::Mat1D m_Weight;

  Matrix::Mat1D m_P0Vec;
  Matrix::Mat1D m_P1Vec;
  float m_Log1PAusuive = 0;
  float m_Log0PAusuive = 0;

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
  void handleSplitPredictData(const char *thBuffer, int startline, int endline);
  void handleSplitTrainData(const char *thBuffer, int pid, long long start,
                            long long end);
  void doTrain();
  inline float sigmod(float x);
};
/***********************************************************
***************************文件处理**************************
************************************************************/
inline float Logistics::sigmod(float x) { return 1.0 / (1 + std::exp(-x)); }
/*
void Logistics::doTrain() {
  ScopeTime t;

  int batch = 8, items = m_samples / batch;
  float rate = 0.02;
  float beta1 = 0.9;
  float beta2 = 0.999;
  float epsilon = 1e-8;

  m_Weight = Matrix::Mat1D(m_features, 0.25);
  Matrix::Mat1D Vdw(m_features, 0);
  Matrix::Mat1D Sdw(m_features, 0);

  for (int epoch = 0; epoch < 10; ++epoch) {
    int start = 0;
    for (int item = 0; item < items; ++item, start += batch) {
      int end = start + batch;

      Matrix::Mat1D grads(m_features, 0);
      for (int i = start; i < end; ++i) {
        float sgd = Dot(m_TrainData[i], m_Weight);
        sgd = sigmod(sgd) - m_Label[i];
        for (int j = 0; j < m_features; ++j) {
          grads[j] += sgd * m_TrainData[i][j];
        }
      }
      for (int i = 0; i < m_features; ++i) {
        float gsv = grads[i] / (float)batch;
        Vdw[i] = beta1 * Vdw[i] + (1 - beta1) * gsv;
        Sdw[i] = beta2 * Sdw[i] + (1 - beta2) * gsv * gsv;
        float cat_vdw = Vdw[i] / (1 - std::pow(beta1, epoch + 1));
        float cat_sdw = Sdw[i] / (1 - std::pow(beta2, epoch + 1));
        m_Weight[i] -= rate * (cat_vdw / std::sqrt(cat_sdw) + epsilon);
      }
    }
  }
  // std::sort(m_Weight.begin(), m_Weight.end());
  std::cerr << "@ train: ";
  t.LogTime();
}
*/
void Logistics::doTrain() {
  ScopeTime t;

  float p0total = 1.0, p1total = 1.0;
  int positive = 0;
  m_P0Vec = Matrix::Mat1D(m_features, 1.0);
  m_P1Vec = Matrix::Mat1D(m_features, 1.0);
  int total = 0;

  for (int pid = 0; pid < NTHREAD; ++pid) {
    const auto &train = m_ThreadData[pid].train;
    const auto &label = m_ThreadData[pid].label;
    const auto &sum = m_ThreadData[pid].sum;
    int up = m_ThreadData[pid].count;
    for (int i = 0; i < up; ++i) {
      if (label[i] == 1) {
        p1total += sum[i];
        Add(m_P1Vec, train[i]);
        ++positive;
      } else {
        p0total += sum[i];
        Add(m_P0Vec, train[i]);
      }
      ++total;
      if (total >= m_samples) {
        break;
      }
    }
    if (total >= m_samples) {
      break;
    }
  }
  float pausuive = (float)(positive) / (float)m_samples;
  m_Log1PAusuive = std::log(pausuive);
  m_Log0PAusuive = std::log(1.0 - pausuive);
  for (auto &it : m_P0Vec) it = std::log(it / p0total);
  for (auto &it : m_P1Vec) it = std::log(it / p1total);

  std::cerr << "@ train: ";
  t.LogTime();
}

void Logistics::handleSplitTrainData(const char *buffer, int pid,
                                     long long start, long long end) {
  Matrix::Mat1D features(m_features);
  int x1, x2, x3, x4;
  int pidx = 0, idx = 0;
  float num, sum = 0;

  int delta = (1000 - m_features) * 6;

  while (start < end) {
    idx = 0;
    sum = 0;
    int j = 0;
    while (true) {
      const char *ptr = &buffer[start + j];
      int k = 0;
      while (k < 60) {
        int x = 1;
        if (ptr[k] == '-') {
          ++k;
          x = -1;
        }
        x1 = ptr[k] - '0';
        x2 = ptr[k + 2] - '0';
        x3 = ptr[k + 3] - '0';
        x4 = ptr[k + 4] - '0';
        num = m_PW[x1][0] + m_PW[x2][1] + m_PW[x3][2] + m_PW[x4][3];
        num *= x;
        features[idx++] = num;
        sum += num;
        k += 6;
        if (idx == m_features) {
          break;
        }
      }
      if (idx == m_features) {
        k += delta;
        while (ptr[k + 1] != '\n') ++k;
        int label = ptr[k] - '0';
        k += 2;
        m_ThreadData[pid].train[pidx] = features;
        m_ThreadData[pid].label[pidx] = label;
        m_ThreadData[pid].sum[pidx] = sum;
        sum = 0;
        idx = 0;
        ++pidx;
        j += k;
        break;
      }
      j += k;
    }
    start += j;
  }
  m_ThreadData[pid].count = pidx;
}

void Logistics::Train() {
  ScopeTime t;
  m_samples = TRAIN_NUM;

  int fd = open(m_trainFile.c_str(), O_RDONLY);
  long long bufsize = m_samples * 7000;
  char *buffer = (char *)mmap(NULL, bufsize, PROT_READ, MAP_PRIVATE, fd, 0);

  int maxItem = m_samples / NTHREAD + 200;
  std::vector<std::thread> Threads(NTHREAD);

  for (int i = 0; i < NTHREAD; ++i) {
    m_ThreadData.emplace_back(Node(maxItem, m_features));
  }

  long long block = bufsize / NTHREAD;
  long long start = 0;
  for (int i = 0; i < NTHREAD; ++i) {
    long long end = std::min(bufsize - 1, start + block);
    if (i == NTHREAD - 1) {
      while (buffer[end] != '\n') --end;
    } else {
      while (buffer[end] != '\n') ++end;
      ++end;
    }
    Threads[i] = std::thread(&Logistics::handleSplitTrainData, this, buffer, i,
                             start, end);
    start = end;
  }
  for (auto &it : Threads) it.join();
  std::cerr << "* TrainData: (" << m_samples << ", " << m_features << ")\n";
  std::cerr << "@ load train: ";
  t.LogTime();

  this->doTrain();
}

void Logistics::handleSplitPredictData(const char *thBuffer, int startline,
                                       int endline) {
  int linesize = 6000;
  long long move = startline * linesize;
  Matrix::Mat1D features(m_features);
  int x1, x2, x3, x4;
  int idx;
  float num;
  for (int i = startline; i < endline; ++i, move += linesize) {
    idx = 0;
    for (int j = 0; j < linesize; j += 60) {
      const char *ptr = &thBuffer[move + j];
      for (int k = 0; k < 60; k += 6) {
        x1 = ptr[k] - '0';
        x2 = ptr[k + 2] - '0';
        x3 = ptr[k + 3] - '0';
        x4 = ptr[k + 4] - '0';
        num = m_PW[x1][0] + m_PW[x2][1] + m_PW[x3][2] + m_PW[x4][3];
        features[idx++] = num;
        if (idx == m_features) {
          break;
        }
      }
      if (idx == m_features) {
        break;
      }
    }
    m_Answer[i] = this->doPredict(features);
  }
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

  // 创建线程
  int start = 0, block = linenum / NTHREAD;
  std::vector<std::thread> Threads(NTHREAD);
  for (int i = 0; i < NTHREAD; ++i) {
    int end = (i == NTHREAD - 1 ? linenum : start + block);
    Threads[i] = std::thread(&Logistics::handleSplitPredictData, this, buffer,
                             start, end);
    start += block;
  }
  for (auto &it : Threads) it.join();

  std::cerr << "@ predict: ";
  t.LogTime();
}

/***********************************************************
***************************训练相关**************************
************************************************************/

void Logistics::Add(Matrix::Mat1D &mat1, const Matrix::Mat1D &mat2) {
  for (int i = 0; i < m_features; i += 8) {
    float *ptr1 = &mat1[i];
    const float *ptr2 = &mat2[i];
    for (int k = 0; k < 8; ++k) {
      ptr1[k] += ptr2[k];
    }
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

int Logistics::doPredict(const Matrix::Mat1D &data) {
  float p1 = Dot(data, m_P1Vec) + m_Log1PAusuive;
  float p0 = Dot(data, m_P0Vec) + m_Log0PAusuive;
  return (p1 > p0 ? 1 : 0);
  // float value = sigmod(Dot(data, m_Weight));
  // return (value > 0.5 ? 1 : 0);
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
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(nullptr);

  const std::string TRAIN = "/data/train_data.txt";
  const std::string PREDICT = "/data/test_data.txt";
  const std::string RESULT = "/projects/student/result.txt";
  const std::string ANSWER = "/projects/student/answer.txt";
  const std::string LOCAL_TRAIN = "/dev/shm/data/train_data.txt";
  const std::string LOCAL_PREDICT = "/dev/shm/data/test_data.txt";
  const std::string LOCAL_RESULT = "/dev/shm/data/result.txt";
  const std::string LOCAL_ANSWER = "/dev/shm/data/answer.txt";

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
