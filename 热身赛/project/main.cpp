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
  typedef std::vector<int> Mat1D;
  typedef std::vector<std::vector<int>> Mat2D;
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

 private:
  int doPredict(const Matrix::Mat1D &data);

 private:
  const int m_samples = 6000;
  const int m_features = 1000;
  const int NTHREAD = 4;  // 线程个数

 private:
  std::vector<Matrix::Mat1D> m_MeanLabel;
  Matrix::Mat1D m_CountLabel;
  Matrix::Mat1D m_Answer;

 private:
  std::string m_trainFile;
  std::string m_predictFile;
  std::string m_resultFile;
  std::string m_answerFile;
};

void Logistics::InitData() {
  m_MeanLabel = Matrix::Mat2D(2, Matrix::Mat1D(m_features, 0));
  m_CountLabel = Matrix::Mat1D(2, 0);
  m_Answer = Matrix::Mat1D(20000, 0);
}
void Logistics::Train() {
  ScopeTime t;

  std::vector<Matrix::Mat2D> threadSum(
      NTHREAD, Matrix::Mat2D(2, Matrix::Mat1D(m_features, 0)));
  Matrix::Mat2D threadCount(NTHREAD, Matrix::Mat1D(2, 0));

  auto foo = [&](int pid, long long start, long long end) {
    std::ifstream thfin(m_trainFile);
    thfin.seekg(start);
    std::string line;

    while (start < end && std::getline(thfin, line)) {
      int label = line.back() - '0';
      ++threadCount[pid][label];
      int pos = 0;
      for (int i = 0; i < m_features; ++i, pos += 6) {
        bool sign = false;
        if (line[pos] == '-') {
          ++pos;
          sign = true;
        }
        int num = (line[pos + 2] - '0') * 100 + (line[pos + 3] - '0') * 10;
        if (sign) num = -num;
        threadSum[pid][label][i] += num;
      }
      start += line.size() + 1;
    }
  };

  std::ifstream fin(m_trainFile);
  long long bufsize = m_samples * 7000;
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
    for (int j = 0; j < m_features; ++j) {
      m_MeanLabel[0][j] += threadSum[i][0][j];
      m_MeanLabel[1][j] += threadSum[i][1][j];
    }
    m_CountLabel[0] += threadCount[i][0];
    m_CountLabel[1] += threadCount[i][1];
  }
  int totalCount = m_CountLabel[0] + m_CountLabel[1];
  for (; totalCount < m_samples; ++totalCount) {
    std::string line;
    std::getline(fin, line);
    int label = line.back() - '0';
    ++m_CountLabel[label];
    int pos = 0;
    for (int i = 0; i < m_features; ++i, pos += 6) {
      bool sign = false;
      if (line[pos] == '-') {
        ++pos;
        sign = true;
      }
      int num = (line[pos + 2] - '0') * 100 + (line[pos + 3] - '0') * 10;
      if (sign) num = -num;
      m_MeanLabel[label][i] += num;
    }
  }
  for (int i = 0; i < m_features; ++i) {
    m_MeanLabel[0][i] /= m_CountLabel[0];
    m_MeanLabel[1][i] /= m_CountLabel[1];
  }

  fin.close();
  std::cerr << "@ train: ";
  t.LogTime();
}
int Logistics::doPredict(const Matrix::Mat1D &data) {
  int dis = 0;
  for (int i = 0; i < m_features; ++i) {
    dis += (data[i] - m_MeanLabel[1][i]) * (data[i] - m_MeanLabel[1][i]) -
           (data[i] - m_MeanLabel[0][i]) * (data[i] - m_MeanLabel[0][i]);
  }
  return dis < 0 ? 1 : 0;
}
void Logistics::Predict() {
  ScopeTime t;

  auto foo = [&](int start, int end) {
    std::ifstream thfin(m_predictFile);
    thfin.seekg((long long)start * 6000L, std::ios::beg);
    std::string line;
    Matrix::Mat1D features(m_features);
    while (start < end && std::getline(thfin, line)) {
      int pos = 0;
      for (int i = 0; i < m_features; ++i, pos += 6) {
        int num = (line[pos + 2] - '0') * 100 + (line[pos + 3] - '0') * 10;
        features[i] = num;
      }
      m_Answer[start++] = doPredict(features);
    }
  };

  std::ifstream fin(m_predictFile);
  long long bufsize = fin.seekg(0, std::ios::end).tellg();
  int linecount = bufsize / 6000;
  int block = linecount / NTHREAD;
  int start = 0, end = block;
  std::vector<std::thread> Threads(NTHREAD);
  for (int i = 0; i < NTHREAD; ++i) {
    end = (i == NTHREAD - 1 ? linecount : start + block);
    Threads[i] = std::thread(foo, start, end);
    start += block;
  }
  for (auto &it : Threads) it.join();

  fin.close();
  std::cerr << "@ predict: ";
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
  const std::string LOCAL_TRAIN = "../data/train_data.txt";
  const std::string LOCAL_PREDICT = "../data/test_data.txt";
  const std::string LOCAL_RESULT = "../data/result.txt";
  const std::string LOCAL_ANSWER = "../data/answer.txt";

  ScopeTime t;

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
  std::cerr << "@ total time: ";
  t.LogTime();
  return 0;
}
