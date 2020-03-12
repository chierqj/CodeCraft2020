#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <thread>
#include <vector>

const std::string TRAIN_FILE = "/data/train_data.txt";
const std::string TEST_FILE = "/data/test_data.txt";
const std::string PREDICT_FILE = "/projects/student/result.txt";
const std::string ANSWER_FILE = "/projects/student/answer.txt";

const std::string LOCAL_TRAIN_FILE = "../data/train_data.txt";
const std::string LOCAL_TEST_FILE = "../data/test_data.txt";
const std::string LOCAL_PREDICT_FILE = "../data/result.txt";
const std::string LOCAL_ANSWER_FILE = "../data/answer.txt";

/***************************************************
 * 打印程序运行时间
 * 1. 函数入口定义变量: ScopeTime t;
 * 2. 函数出口打印时间: t.LogTime();
 **************************************************/
class ScopeTime {
 public:
  ScopeTime() : m_begin(std::chrono::high_resolution_clock::now()) {}
  void LogTime() const {
    auto t = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - m_begin);
    double elapsed = (double)(t.count() * 1.0) / 1000.0;
    std::cerr << "* 运行时间: " << elapsed << "s\n";
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

/***************************************************
 * 矩阵运算
 * 1. NewMat1D(n, value)生成一数组，默认value
 * 2. NewMat2D(n, m, value)生成二维数组，默认value
 * 3. T(Mat2D mat) mat.T
 * 4. Subtraction(Mat1D mat1, Mat1D mat2) mat1 - mat2
 * 5. Dot(Mat1D mat1, Mat1D mat2) mat1(1*n) 点乘 mat2(1*n)
 * 6. Dot(Mat2D mat1, Mat1D mat2) mat1(n*m) 点乘 mat2(1*m)
 * 7. Multipy(Mat2D mat1, Mat2D mat2) mat1(n*m)*mat2(m*k)
 **************************************************/
struct Matrix {
  typedef std::vector<double> Mat1D;
  typedef std::vector<std::vector<double>> Mat2D;
};
// 一维 + 一维
Matrix::Mat1D operator+(const Matrix::Mat1D &mat1, const Matrix::Mat1D &mat2) {
  int n1 = mat1.size(), n2 = mat2.size();
  assert(n1 == n2);
  Matrix::Mat1D ans = Matrix::Mat1D(n1, 0);
  for (int i = 0; i < n1; ++i) {
    ans[i] = mat1[i] + mat2[i];
  }
  return ans;
}
// 二维 + 一维
Matrix::Mat2D operator+(const Matrix::Mat2D &mat1, const Matrix::Mat1D &mat2) {
  Matrix::Mat2D ans;
  for (auto &it : mat1) {
    ans.emplace_back(it + mat2);
  }
  return ans;
}
// 一维 - 一维
Matrix::Mat1D operator-(const Matrix::Mat1D &mat1, const Matrix::Mat1D &mat2) {
  int n1 = mat1.size(), n2 = mat2.size();
  assert(n1 == n2);
  Matrix::Mat1D ans = Matrix::Mat1D(n1, 0);
  for (int i = 0; i < n1; ++i) {
    ans[i] = mat1[i] - mat2[i];
  }
  return ans;
}
// 一维 - 二维
Matrix::Mat2D operator-(const Matrix::Mat1D &mat1, const Matrix::Mat2D &mat2) {
  Matrix::Mat2D ans;
  for (auto &it : mat2) {
    ans.emplace_back(mat1 - it);
  }
  return ans;
}
// 二维 - 一维
Matrix::Mat2D operator-(const Matrix::Mat2D &mat1, const Matrix::Mat1D &mat2) {
  Matrix::Mat2D ans;
  for (auto &it : mat1) {
    ans.emplace_back(it - mat2);
  }
  return ans;
}
// 转置
Matrix::Mat2D T(const Matrix::Mat2D &mat) {
  int n = mat.size();
  assert(n > 0);
  int m = mat[0].size();
  assert(m > 0);
  Matrix::Mat2D ans(m, Matrix::Mat1D(n));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      ans[j][i] = mat[i][j];
    }
  }
  return ans;
}
// 一维 点乘 一维
double Dot(const Matrix::Mat1D &mat1, const Matrix::Mat1D &mat2) {
  int n1 = mat1.size(), n2 = mat2.size();
  assert(n1 == n2);
  double ans = 0.0L;
  for (int i = 0; i < n1; ++i) {
    ans += mat1[i] * mat2[i];
  }
  return ans;
}
// 二维 点乘 一维
Matrix::Mat1D Dot(const Matrix::Mat2D &mat1, const Matrix::Mat1D &mat2) {
  Matrix::Mat1D ans;
  for (auto &it : mat1) {
    ans.emplace_back(Dot(it, mat2));
  }
  return ans;
}
// 一 点乘 二维
Matrix::Mat1D Dot(const Matrix::Mat1D &mat1, const Matrix::Mat2D &mat2) {
  int n1 = mat1.size(), n2 = mat2.size();
  assert(n2 > 0);
  int m2 = mat2[0].size();
  assert(n1 == n2);
  Matrix::Mat1D ans(m2, 0.0L);
  for (int j = 0; j < m2; ++j) {
    for (int i = 0; i < n2; ++i) {
      ans[j] += mat1[i] * mat2[i][j];
    }
  }
  return ans;
}
static Matrix::Mat2D operator*(const Matrix::Mat2D &mat1,
                               const Matrix::Mat2D &mat2) {
  Matrix::Mat2D mat2T = T(mat2);
  int n = mat1.size(), m = mat2[0].size();
  Matrix::Mat2D mat(n, Matrix::Mat1D(m));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      mat[i][j] = Dot(mat1[i], mat2T[j]);
    }
  }
  return mat;
};
static Matrix::Mat1D operator*(const Matrix::Mat2D &mat1,
                               const Matrix::Mat1D &mat2) {
  int n = mat1.size();
  Matrix::Mat1D mat;
  for (const auto &x : mat1) {
    mat.emplace_back(Dot(x, mat2));
  }
  return mat;
};

/***************************************************
 * Model
 * 1. 构造函数Model(Matrix::Mat2D, Matrix::Mat1D) 传入train和label
 * 2. Train() 训练模型
 * 3. Predict(Matrix::Mat1D) 预测数据必须经过预处理，特征新加一列1
 * 4. Print() 打印相关系数
 * 5. sigmod(double x) 计算激活函数
 * 6. splitTrainData() 按照设置的比率拆分训练集，训练集+验证集
 * 7. score() 对验证集进行评分
 **************************************************/
class Model {
 public:
  static const int MAX_ITER_TIME = 120;     // 迭代次数
  const int SKIP_SAMPLES = 2;               // 随机梯度下降(样本x选1)
  static const int SELECT_TRAIN_NUM = 800;  // 选择样本数(-1表示全选)
  const double WEIGHT = 1.0;                // 初始化weight
  const int SHOW_TRAIN_STEP = 50;           // 每隔多少代打印log
  const double PREDICT_TRUE_THRESH = 0.5;   // 划分答案

 public:
  void Train(const Matrix::Mat2D &Train, const Matrix::Mat1D &Label);
  inline int Predict(const Matrix::Mat1D &data);

 private:
  inline double sigmod(const double &x);  // sigmod
  void gd(const Matrix::Mat2D &Train, const Matrix::Mat1D &Label,
          const int &epoch);
  void sgd(const Matrix::Mat2D &Train, const Matrix::Mat1D &Label,
           std::vector<int> &RandomIndex, const int &epoch);

 private:
  Matrix::Mat1D m_Weight;  // 参数
  int m_samples = 0;
  int m_features = 0;
};
inline int Model::Predict(const Matrix::Mat1D &data) {
  double sigValue = this->sigmod(Dot(data, m_Weight));
  return (sigValue >= PREDICT_TRUE_THRESH ? 1 : 0);
}
inline double Model::sigmod(const double &x) {
  return 1.0 / (1.0 + std::exp(-x));
}
void Model::gd(const Matrix::Mat2D &Train, const Matrix::Mat1D &Label,
               const int &epoch) {
  static Matrix::Mat2D trainT = T(Train);
  double learningRate = 5.0 / (epoch + 1) + 0.001;
  auto mat = Train * m_Weight;
  for (auto &x : mat) {
    x = sigmod(x);
  }
  auto sig = mat - Label;
  sig = trainT * sig;
  for (int i = 0; i < m_features; ++i) {
    m_Weight[i] -= sig[i] * learningRate;
  }
}
void Model::sgd(const Matrix::Mat2D &Train, const Matrix::Mat1D &Label,
                std::vector<int> &RandomIndex, const int &epoch) {
  Tools::ShuffleVector(RandomIndex);
  for (int j = 0; j < m_samples; ++j) {
    if (j % SKIP_SAMPLES == 0) {  // SKIP_SAMPLES选一个
      int index = RandomIndex[j];
      double alpha = 4.0 / (epoch + j + 1) + 0.01;
      double sgd = this->sigmod(Dot(Train[index], m_Weight));
      sgd -= Label[index];
      for (int k = 0; k < m_features; ++k) {
        m_Weight[k] -= alpha * sgd * Train[index][k];
      }
    }
  }
}
void Model::Train(const Matrix::Mat2D &Train, const Matrix::Mat1D &Label) {
  std::cerr << "--------------------------------------\n";
  std::cerr << "* 开始训练\n";
  ScopeTime t;

  m_samples = Train.size();
  m_features = Train[0].size();

  m_Weight = Matrix::Mat1D(m_features, WEIGHT);
  m_Weight.back() = 0.0;

  std::vector<int> RandomIndex;
  for (int i = 0; i < Train.size(); ++i) {
    RandomIndex.emplace_back(i);
  }

  for (int epoch = 0; epoch < MAX_ITER_TIME; ++epoch) {
    // this->sgd(Train, Label, RandomIndex, epoch);
    this->gd(Train, Label, epoch);

    if (epoch % SHOW_TRAIN_STEP == 0) {
      std::cerr << "* epoch: " << epoch << "\n";
    }
  }
  t.LogTime();
  std::cerr << "--------------------------------------\n";
}

/***************************************************
 * Adama优化器
 **************************************************/
class Adam {
 public:
  const int MAX_ITER_TIME = 1000;

 public:
  void Train(const Matrix::Mat2D &Train, const Matrix::Mat1D &Label);
  double Predict(const Matrix::Mat1D &data);

 private:
  inline double sigmod(double z);
  inline Matrix::Mat1D sigmod(const Matrix::Mat1D &mat);

 private:
  Matrix::Mat1D m_Weight;
  Matrix::Mat1D m_B;

  double m_m = 0;
  double m_b = 0;
  double m_alpha = 0.1;
  double m_threshold = 0.5;
  double m_batch = 128;
  double m_beta1 = 0.9;
  double m_beta2 = 0.999;
  double m_epsilon = 1e-8;
};
inline double Adam::sigmod(double z) { return 1.0 / (1.0 + std::exp(-z)); }
inline Matrix::Mat1D Adam::sigmod(const Matrix::Mat1D &mat) {
  Matrix::Mat1D ans;
  for (auto &v : mat) {
    ans.emplace_back(this->sigmod(v));
  }
  return ans;
}
void Adam::Train(const Matrix::Mat2D &Train, const Matrix::Mat1D &Label) {
  int n_samples = Train.size(), n_features = Train[0].size();
  auto Vdw = Matrix::Mat1D(n_features, 0);
  double Vdb = 0;
  auto Sdw = Matrix::Mat1D(n_features, 0);
  double Sdb = 0;
  int batchNum = m_m / m_batch + 1;
}

/***************************************************
 * Simulation
 * 1. LoadData() 加载训练集和预测集
 * 2. LoadAnswer() 如果是本地的话，加载answer
 * 3. Train() 数据预处理+模型训练
 * 4. Predict() 根据m_model提供的Predict接口预测数据
 * 5. Score() 根据本地answer评分
 * 6. SaveAnswer() 保存答案
 * 7. calculateGini() 计算每一维度的gini系数
 **************************************************/
class Simulation {
 private:
  Model *m_model;

 public:
  Simulation();
  ~Simulation();
  void LoadData();
  void LoadAnswer();
  void Train();
  void PredictAndSaveAnswer();
  void Score();
  void SaveAnswer();

 private:
  void loadDataFromCharVec(const std::vector<char> &chars, Matrix::Mat2D &ans,
                           bool trainfile);
  void loadTrainByCharVec();
  void loadTrainBySkipNumber(int skip);
  void loadTestByCharVec();

 private:
  std::string m_trainFile;        // 训练集路径
  std::string m_testFile;         // 预测集路径
  std::string m_predictFile;      // 结果保留路径
  std::string m_answerFile;       // 正确答案路径
  Matrix::Mat2D m_TrainData;      // 训练集合
  Matrix::Mat2D m_PredictData;    // 预测集合
  std::vector<int> m_AnswerData;  // 正确答案
  int m_positiveSamples = 0;      // 正样本个数
  int m_negativeSamples = 0;      // 负样本个数
  std::string m_answer;           // 答案
};

Simulation::Simulation() {
#ifdef LOCAL_TRAIN
  m_trainFile = LOCAL_TRAIN_FILE;
  m_testFile = LOCAL_TEST_FILE;
  m_predictFile = LOCAL_PREDICT_FILE;
  m_answerFile = LOCAL_ANSWER_FILE;
#else
  m_trainFile = TRAIN_FILE;
  m_testFile = TEST_FILE;
  m_predictFile = PREDICT_FILE;
  m_answerFile = ANSWER_FILE;
#endif
}
Simulation::~Simulation() { delete m_model; }

void Simulation::loadDataFromCharVec(const std::vector<char> &chars,
                                     Matrix::Mat2D &ans, bool trainfile) {
  bool negative = false, dot = false;
  double minx = 0, maxx = 1.0;
  double integer = 0, decimal = 0, multiple = 10;
  Matrix::Mat1D features;
  for (auto &v : chars) {
    if (v == ',' || v == '\n') {
      double num = integer + decimal;
      if (negative) num = -num;
      if (num < 0) num = 0;
      if (num > 1) num = 1;
      negative = dot = false;
      integer = decimal = 0;
      multiple = 10;

      features.emplace_back(num);

      if (v == ',') {
        minx = std::min(minx, num);
        maxx = std::max(maxx, num);
      } else {
        if (trainfile) {
          if (minx >= 0 && maxx <= 1) {
            ans.emplace_back(features);
          }
        } else {
          ans.emplace_back(features);
        }
        features.clear();
        minx = 0;
        maxx = 1.0;
      }
    } else if (v == '-') {
      negative = true;
    } else if (v == '.') {
      dot = true;
    } else {
      double x = v - '0';
      if (!dot) {
        integer = integer * 10 + x;
      } else {
        decimal = decimal + x / multiple;
        multiple *= 10;
      }
    }
  }
}

void Simulation::loadTrainByCharVec() {
  std::ifstream fin(m_trainFile, std::ios::binary);
  assert(fin);
  std::vector<char> buf(
      static_cast<unsigned int>(fin.seekg(0, std::ios::end).tellg()));
  fin.seekg(0, std::ios::beg)
      .read(&buf[0], static_cast<std::streamsize>(buf.size()));
  fin.close();
  this->loadDataFromCharVec(buf, m_TrainData, true);
}

void Simulation::loadTestByCharVec() {
  std::ifstream fin(m_testFile, std::ios::binary);
  assert(fin);
  std::vector<char> buf(
      static_cast<unsigned int>(fin.seekg(0, std::ios::end).tellg()));
  fin.seekg(0, std::ios::beg)
      .read(&buf[0], static_cast<std::streamsize>(buf.size()));
  fin.close();
  this->loadDataFromCharVec(buf, m_PredictData, false);
}

double ToDouble(std::string &s) {
  int n = s.size(), i = n - 1;
  double ans = 0;
  if (s[0] == '-') {
    return 0;
  } else if (s[0] == '0') {
    while (i > 1) {
      ans = ans + (s[i--] - '0');
      ans /= 10;
    }
    return ans;
  } else {
    return 1;
  }
}
void Simulation::loadTrainBySkipNumber(int skip) {
  std::ifstream fin(m_trainFile);
  std::string line, temp;
  int cnt = 0;
  while (fin) {
    if (cnt >= skip) break;
    std::getline(fin, line);
    if (line.empty()) {
      continue;
    }
    std::istringstream iss(line);
    std::vector<double> features;
    bool sign = true;
    while (getline(iss, temp, ',')) {
      double num = ToDouble(temp);
      features.emplace_back(num);
    }
    if (sign) {
      m_TrainData.emplace_back(features);
      ++cnt;
    }
  }
  fin.close();
}

void Simulation::LoadAnswer() {
#ifndef LOCAL_TRAIN
  return;
#endif
  std::ifstream fin(m_answerFile);
  assert(fin);
  int x;
  while (fin >> x) {
    m_AnswerData.emplace_back(x);
  }
  fin.close();
}

void Simulation::Train() {
  Matrix::Mat1D label;
  for (auto &it : m_TrainData) {
    int n_label = it.back();
    it.back() = 1;
    if (n_label == 1) {
      ++m_positiveSamples;
    } else {
      ++m_negativeSamples;
    }
    label.emplace_back(n_label);
  }
  m_model = new Model();
  m_model->Train(m_TrainData, label);
}

void Simulation::PredictAndSaveAnswer() {
  for (auto &test : m_PredictData) {
    test.emplace_back(1.0);
    char c = m_model->Predict(test) + '0';
    m_answer.push_back(c);
    m_answer.push_back('\n');
  }
  std::ofstream fout(m_predictFile);
  auto &rdbuf(*fout.rdbuf());
  rdbuf.sputn(m_answer.data(), m_answer.size());
  fout.close();
}

void Simulation::Score() {
#ifndef LOCAL_TRAIN
  return;
#endif
  std::cerr << "--------------------------------------\n";
  std::cerr << "* 开始评分\n";

  int sz = m_PredictData.size();
  int same = 0, unsame = 0;
  for (int i = 0; i < sz; ++i) {
    if (m_answer[i * 2] - '0' == m_AnswerData[i]) {
      ++same;
    } else {
      ++unsame;
    }
  }
  // std::cerr << "* 正样本: " << m_positiveSamples << "\n";
  // std::cerr << "* 负样本: " << m_negativeSamples << "\n";
  // std::cerr << "* 迭代数: " << Model::MAX_ITER_TIME << "\n";
  // std::cerr << "* 正确数: " << same << "\n";
  // std::cerr << "* 错误数: " << unsame << "\n";
  std::cerr << "* 正确率: " << (double)same / (double)sz * 100 << "%\n";
  std::cerr << "--------------------------------------\n";
}

void Simulation::LoadData() {
  std::cerr << "--------------------------------------\n";
  std::cerr << "* 加载数据\n";
  ScopeTime t;
  this->loadTrainBySkipNumber(Model::SELECT_TRAIN_NUM);
  this->loadTestByCharVec();
  this->LoadAnswer();
  std::cerr << "* TrainData: (" << m_TrainData.size() << ", "
            << m_TrainData[0].size() << ")\n";
  std::cerr << "* TestData: (" << m_PredictData.size() << ", "
            << m_PredictData[0].size() << ")\n";
  t.LogTime();
  std::cerr << "--------------------------------------\n";
}

int main() {
  std::cerr << std::fixed << std::setprecision(3);

  ScopeTime t;
  Simulation *simulation = new Simulation();

  simulation->LoadData();
  simulation->Train();
  simulation->PredictAndSaveAnswer();
  simulation->Score();

  t.LogTime();
  delete simulation;
  return 0;
}