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
#include <vector>

const std::string TRAIN_FILE = "/data/train_data.txt";
const std::string TEST_FILE = "/data/test_data.txt";
const std::string PREDICT_FILE = "/projects/student/result.txt";
const std::string ANSWER_FILE = "/projects/student/answer.txt";

const std::string LOCAL_TRAIN_FILE = "./data/train_data.txt";
const std::string LOCAL_TEST_FILE = "./data/test_data.txt";
const std::string LOCAL_PREDICT_FILE = "./data/result.txt";
const std::string LOCAL_ANSWER_FILE = "./data/answer.txt";

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
class Martix {
 public:
  typedef std::vector<double> Mat1D;
  typedef std::vector<std::vector<double>> Mat2D;

  static void Print(const Mat1D &mat) {
    std::cerr << "{";
    for (int i = 0; i < mat.size(); ++i) {
      if (i != 0) {
        std::cerr << ", ";
      }
      std::cerr << mat[i];
    }
    std::cerr << "}\n";
  }
  static void Print(const Mat2D &mat) {
    for (int i = 0; i < mat.size(); ++i) {
      std::cerr << i << ": ";
      Martix::Print(mat[i]);
    }
  }

  static Mat1D NewMat1D(int n, double value);
  static Mat2D NewMat2D(int n, int m, double value);
  static Mat2D T(const Mat2D &mat);
  static Mat1D Subtraction(const Mat1D &mat1, const Mat1D &mat2);

  static double Dot(const Mat1D &mat1, const Mat1D &mat2);
  static Mat1D Dot(const Mat2D &mat2D, const Mat1D &mat1D);
  static Mat2D Multipy(const Mat2D &mat1, const Mat2D &mat2);
};

Martix::Mat1D Martix::NewMat1D(int n, double value) {
  Mat1D mat(n, value);
  return mat;
}
Martix::Mat2D Martix::NewMat2D(int n, int m, double value) {
  Mat2D mat(n, Mat1D(m, value));
  return mat;
}
Martix::Mat2D Martix::T(const Mat2D &mat) {
  int n = mat.size();
  assert(n > 0);
  int m = mat[0].size();
  assert(m > 0);
  Mat2D ans = Martix::NewMat2D(m, n, 0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      ans[j][i] = mat[i][j];
    }
  }
  return ans;
}
Martix::Mat1D Martix::Subtraction(const Mat1D &mat1, const Mat1D &mat2) {
  int n1 = mat1.size(), n2 = mat2.size();
  assert(n1 == n2);
  Mat1D ans = Martix::NewMat1D(n1, 0);
  for (int i = 0; i < n1; ++i) {
    ans[i] = mat1[i] - mat2[i];
  }
  return ans;
}
double Martix::Dot(const Mat1D &mat1, const Mat1D &mat2) {
  int n1 = mat1.size(), n2 = mat2.size();
  assert(n1 == n2);
  double ans = 0.0L;
  for (int i = 0; i < n1; ++i) {
    ans += mat1[i] * mat2[i];
  }
  return ans;
}
Martix::Mat1D Martix::Dot(const Martix::Mat2D &mat2D,
                          const Martix::Mat1D &mat1D) {
  int n = mat2D.size();
  assert(n > 0);
  int m = mat2D[0].size();
  int m1 = mat1D.size();
  assert(m == m1);

  Martix::Mat1D ans;
  for (const auto &row : mat2D) {
    ans.emplace_back(Martix::Dot(row, mat1D));
  }
  return ans;
}
Martix::Mat2D Martix::Multipy(const Martix::Mat2D &mat1,
                              const Martix::Mat2D &mat2) {
  int n1 = mat1.size();
  assert(n1 > 0);
  int m1 = mat1[0].size();
  assert(m1 > 0);
  int n2 = mat2.size();
  assert(n2 > 0);
  int m2 = mat2[0].size();
  assert(m2 > 0);
  assert(m1 == n2);
  auto ans = Martix::NewMat2D(n1, m2, 0);
  for (int i = 0; i < n1; ++i) {
    for (int j = 0; j < m2; ++j) {
      for (int k = 0; k < m1; ++k) {
        ans[i][j] += mat1[i][k] * mat2[k][j];
      }
    }
  }
  return ans;
}

/***************************************************
 * Model
 * 1. 构造函数Model(Martix::Mat2D, Martix::Mat1D) 传入train和label
 * 2. Train() 训练模型
 * 3. Predict(Martix::Mat1D) 预测数据必须经过预处理，特征新加一列1
 * 4. Print() 打印相关系数
 * 5. sigmod(double x) 计算激活函数
 * 6. splitTrainData() 按照设置的比率拆分训练集，训练集+验证集
 * 7. score() 对验证集进行评分
 **************************************************/
class Model {
 public:
  static const int SELECT_FEATURE_NUM = -1;  // 随机选取几个维度(-1全选)
  static const int MAX_ITER_TIME = 100;      // 迭代次数
  static const int SELECT_TRAIN_NUM = 2500;  // 选择样本数(-1表示全选)
  const double TRAIN_SCALE = 1;              // 训练集占比
  const int SKIP_SAMPLES = 1;                // 随机梯度下降(样本x选1)
  const double WEIGHT = 1.0;                 // 初始化weight
  const int SHOW_TRAIN_STEP = 20;            // 每隔多少代打印log
  const double PREDICT_TRUE_THRESH = 0.5;    // 划分答案

 public:
  Model(const Martix::Mat2D &trainSet, const Martix::Mat1D &label);
  void Train();
  int Predict(const Martix::Mat1D &data);
  void Print();

 private:
  inline double sigmod(const double &x);  // sigmod
  void splitTrainData();                  // 拆分训练集合
  double score();                         // 评分

 private:
  Martix::Mat1D m_Weight;     // 参数
  Martix::Mat2D m_TrainData;  // 训练集
  Martix::Mat1D m_Label;      // 标签
  Martix::Mat2D m_TestData;   // 训练集
  Martix::Mat1D m_TestLabel;  // 标签
  int m_samples = 0;          // 数据集大小
  int m_features = 0;         // 特征个数
  int m_testSamples = 0;      // 测试集row
  int m_testFeatures = 0;     // 测试集feature
};

Model::Model(const Martix::Mat2D &trainSet, const Martix::Mat1D &label) {
  m_TrainData = trainSet;
  m_Label = label;
  for (auto &row : m_TrainData) {
    row.emplace_back(1);
  }
  m_samples = m_TrainData.size();
  m_features = m_TrainData[0].size();
}

void Model::Print() {
  std::cerr << "* 训练集: (" << m_samples << ", ";
  std::cerr << m_features << ")\n";
  std::cerr << "* 测试集: (" << m_testSamples << ", ";
  std::cerr << m_testFeatures << ")\n";
}

void Model::Train() {
  this->splitTrainData();
  std::cerr << "--------------------------------------\n";
  std::cerr << "* 开始训练\n";
  ScopeTime t;

  std::vector<int> RandomIndex;
  for (int i = 0; i < m_samples; ++i) {
    RandomIndex.emplace_back(i);
  }
  m_Weight = Martix::NewMat1D(m_features, WEIGHT);
  m_Weight.back() = 0.0;

  for (int i = 0; i < MAX_ITER_TIME; ++i) {
    Tools::ShuffleVector(RandomIndex);
    for (int j = 0; j < m_samples; ++j) {
      if (j % SKIP_SAMPLES == 0) {  // 随机每10个里面选一个进行迭代
        int index = RandomIndex[j];
        double alpha = 4.0 / (i + j + 1) + 0.01;
        double sgd = this->sigmod(Martix::Dot(m_TrainData[index], m_Weight));
        sgd -= m_Label[index];
        for (int k = 0; k < m_features; ++k) {
          m_Weight[k] -= alpha * sgd * m_TrainData[index][k];
        }
      }
    }

#ifdef LOCAL_TRAIN
    double score = this->score();
    if (i % SHOW_TRAIN_STEP == 0) {
      std::cerr << "* 代数" << i << ", 正确率: " << score << "\n";
    }
#endif
  }
  t.LogTime();
  std::cerr << "--------------------------------------\n";
}

int Model::Predict(const Martix::Mat1D &data) {
  double sigValue = this->sigmod(Martix::Dot(data, m_Weight));
  return (sigValue >= PREDICT_TRUE_THRESH ? 1 : 0);
}

inline double Model::sigmod(const double &x) {
  return 1.0 / (1.0 + std::exp(-x));
}

void Model::splitTrainData() {
  if (fabs(TRAIN_SCALE - 1.0) <= 1e-5) {
    return;
  }
  std::cerr << "--------------------------------------\n";
  std::cerr << "* 分割数据\n";
  std::vector<int> zeroVec;
  std::vector<int> oneVec;
  for (int i = 0; i < m_TrainData.size(); ++i) {
    if (m_Label[i] == 0) {
      zeroVec.emplace_back(i);
    } else {
      oneVec.emplace_back(i);
    }
  }

  Martix::Mat2D splitTrainData;
  Martix::Mat1D splitLabel;

  auto foo = [&](std::vector<int> &vt) {
    Tools::ShuffleVector(vt);
    int cutIndex = vt.size() * TRAIN_SCALE;
    for (int i = 0; i < cutIndex; ++i) {
      splitTrainData.emplace_back(m_TrainData[vt[i]]);
      splitLabel.emplace_back(m_Label[vt[i]]);
    }
    for (int i = cutIndex; i < vt.size(); ++i) {
      m_TestData.emplace_back(m_TrainData[vt[i]]);
      m_TestLabel.emplace_back(m_Label[vt[i]]);
    }
  };

  foo(zeroVec);
  foo(oneVec);

  m_TrainData.clear();
  m_Label.clear();
  m_TrainData = splitTrainData;
  m_Label = splitLabel;
  m_samples = m_TrainData.size();
  m_features = m_TrainData[0].size();
  m_testSamples = m_TestData.size();
  if (m_testSamples > 0) {
    m_testFeatures = m_TestData[0].size();
  }

  this->Print();
  std::cerr << "--------------------------------------\n";
}

double Model::score() {
  if (std::fabs(TRAIN_SCALE - 1.0) <= 1e-5) {
    return -1;
  }
  int ac = 0, wa = 0;
  for (int i = 0; i < m_testSamples; ++i) {
    int predictLabel = this->Predict(m_TestData[i]);
    if (predictLabel == m_TestLabel[i]) {
      ++ac;
    } else {
      ++wa;
    }
  }
  double ans = (double)(ac * 1.0) / (double)(ac + wa);
  return ans;
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
struct Data {
  Martix::Mat1D features;
  int label;
};
struct Feature {
  int index;
  std::vector<std::pair<double, int>> values;
  double gain;
};
class Simulation {
 private:
  Model *m_model;

 public:
  Simulation();
  ~Simulation();
  void LoadData();
  void LoadAnswer();
  void Train();
  void Predict();
  void Score();
  void SaveAnswer();

 private:
  void doLoadData(const std::string &file, std::vector<Data> &vt);
  void loadFromCharVec(const std::vector<char> &vt, std::vector<Data> &data,
                       const bool &trainFile);
  void doLoadAnswer();
  void calculateGain();
  Martix::Mat1D getSelectFeatures(const Data &data);

 private:
  std::string m_trainFile;               // 训练集路径
  std::string m_testFile;                // 预测集路径
  std::string m_predictFile;             // 结果保留路径
  std::string m_answerFile;              // 正确答案路径
  std::vector<Data> m_TrainData;         // 训练集合
  std::vector<Data> m_PredictData;       // 预测集合
  std::vector<int> m_AnswerData;         // 正确答案
  std::vector<Feature> m_Features;       // 特征处理
  std::vector<int> m_allowFeatureIndex;  // 允许选择的特征
  int m_positiveSamples = 0;             // 正样本个数
  int m_negativeSamples = 0;             // 负样本个数
  const bool FILTER = true;              // 是否过滤数据[0-1]
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

void Simulation::loadFromCharVec(const std::vector<char> &chars,
                                 std::vector<Data> &data,
                                 const bool &trainFile) {
  bool negative = false;
  bool dot = false;
  double x1 = 0, x2 = 0, r = 10;
  Martix::Mat1D features;
  double minx = 0, maxx = 1.0;
  int cnt = 0;
  for (auto &v : chars) {
    if (v == ',' || v == '\n') {
      double num = x1 + x2;
      if (negative) {
        num = -num;
      }
      negative = false;
      dot = false;
      x1 = 0;
      x2 = 0;
      r = 10;

      if (v == ',') {
        minx = std::min(minx, num);
        maxx = std::max(maxx, num);
        features.emplace_back(num);
      } else {
        if (trainFile) {
          if (minx >= 0 && maxx <= 1.0) {
            int label = num;
            if (label == 1) {
              ++m_positiveSamples;
            } else {
              ++m_negativeSamples;
            }
            data.emplace_back(Data{features, label});
            ++cnt;
            if (Model::SELECT_TRAIN_NUM != -1 &&
                cnt >= Model::SELECT_TRAIN_NUM) {
              return;
            }
          }
        } else {
          features.emplace_back(num);
          data.emplace_back(Data{features, -1});
        }
        features.clear();
        minx = 0;
        maxx = 1.0;
      }
      continue;
    }

    if (v == '-') {
      negative = true;
      continue;
    }
    if (v == '.') {
      dot = true;
      continue;
    }
    double x = v - '0';
    if (!dot) {
      x1 = x1 * 10 + x;
    } else {
      x2 = x2 + x / r;
      r *= 10;
    }
  }
}

void Simulation::doLoadData(const std::string &file, std::vector<Data> &vt) {
  std::ifstream fin(file, std::ios::binary);
  assert(fin);
  std::vector<char> buf(
      static_cast<unsigned int>(fin.seekg(0, std::ios::end).tellg()));
  fin.seekg(0, std::ios::beg)
      .read(&buf[0], static_cast<std::streamsize>(buf.size()));
  fin.close();
  bool sign = (file == m_trainFile);
  this->loadFromCharVec(buf, vt, sign);
}

void Simulation::doLoadAnswer() {
#ifdef LOCAL_TRAIN
  std::ifstream fin(m_answerFile);
  assert(fin);
  int x;
  while (fin >> x) {
    m_AnswerData.emplace_back(x);
  }
  fin.close();
#endif
}

void Simulation::calculateGain() {
  if (Model::SELECT_FEATURE_NUM == -1) {
    return;
  }

  ScopeTime t;
  int total = m_TrainData.size();
  std::cerr << "--------------------------------------\n";
  std::cerr << "* 计算Gain系数\n";
  std::cerr << "* 总样本: " << total << "\n";
  std::cerr << "* 正样本: " << m_positiveSamples << "\n";
  std::cerr << "* 负样本: " << m_negativeSamples << "\n";

  int n = m_TrainData.size(), m = m_TrainData[0].features.size();
  for (int i = 0; i < m; ++i) {
    std::vector<std::pair<double, int>> values;
    for (auto &train : m_TrainData) {
      values.emplace_back(std::make_pair(train.features[i], train.label));
    }
    std::sort(values.begin(), values.end());
    m_Features.emplace_back(Feature{i, values, 0.0});
  }

  auto getEnt = [](double x, double y) {
    double e1 = (x == 0 ? 0 : x * std::log2(x));
    double e2 = (y == 0 ? 0 : y * std::log2(y));
    return -(e1 + e2);
  };

  double posRate = (double)(m_positiveSamples) / total;
  double negRate = 1.0 - posRate;
  double Ent = getEnt(posRate, negRate);

  for (auto &fea : m_Features) {
    std::vector<int> prefix(total, 0);
    int cnt = 0;
    double maxGain = 0;
    for (int i = 0; i < total; ++i) {
      if (fea.values[i].second == 1) {
        ++cnt;
      }
      prefix[i] = cnt;
    }
    for (int i = 1; i < total; ++i) {
      double x1 = (double)prefix[i - 1] / i;
      double y1 = 1.0 - x1;
      double x2 = (double)(prefix[total - 1] - prefix[i - 1]) / (total - i);
      double y2 = 1.0 - x2;
      double ent1 = getEnt(x1, y1), ent2 = getEnt(x2, y2);
      double rate = (double)i / total;
      double gain = Ent - (rate * ent1 + (1 - rate) * ent2);
      maxGain = std::max(maxGain, gain);
    }
    fea.gain = maxGain;
  }
  std::sort(m_Features.begin(), m_Features.end(),
            [&](const Feature &fea1, const Feature &fea2) {
              return fea1.gain > fea2.gain;
            });

  t.LogTime();
  std::cerr << "--------------------------------------\n";
}

Martix::Mat1D Simulation::getSelectFeatures(const Data &data) {
  static int sz = m_TrainData[0].features.size();
  Martix::Mat1D fea;
  int up = std::min(Model::SELECT_FEATURE_NUM, sz);
  for (int i = 0; i < up; ++i) {
    int index = m_Features[i].index;
    fea.emplace_back(m_TrainData[i].features[index]);
  }
  return fea;
}

void Simulation::Train() {
  this->calculateGain();

  Martix::Mat2D train;
  Martix::Mat1D label;

  if (Model::SELECT_FEATURE_NUM != -1) {
    for (auto &it : m_TrainData) {
      Martix::Mat1D fea = this->getSelectFeatures(it);
      train.emplace_back(fea);
      label.emplace_back(it.label);
    }
  } else {
    for (auto &it : m_TrainData) {
      train.emplace_back(it.features);
      label.emplace_back(it.label);
    }
  }

  m_model = new Model(train, label);
  m_model->Train();
}

void Simulation::Predict() {
  for (auto &test : m_PredictData) {
    Martix::Mat1D fea;
    if (Model::SELECT_FEATURE_NUM != -1) {
      fea = this->getSelectFeatures(test);
    } else {
      fea = test.features;
    }
    fea.emplace_back(1.0);
    test.label = m_model->Predict(fea);
  }
}

void Simulation::Score() {
#ifdef LOCAL_TRAIN
  std::cerr << "--------------------------------------\n";
  std::cerr << "* 开始评分\n";

  int sz = m_PredictData.size();
  int same = 0, unsame = 0;
  for (int i = 0; i < sz; ++i) {
    if (m_PredictData[i].label == m_AnswerData[i]) {
      ++same;
    } else {
      ++unsame;
    }
  }
  m_model->Print();
  std::cerr << "* 迭代: " << Model::MAX_ITER_TIME << "\n";
  std::cerr << "* 正确: " << same << "\n";
  std::cerr << "* 错误: " << unsame << "\n";
  std::cerr << "* 正确率: " << (double)same / (double)sz << "\n";
  std::cerr << "--------------------------------------\n";
#endif
}

void Simulation::LoadData() {
  std::cerr << "--------------------------------------\n";
  std::cerr << "* 加载数据\n";
  ScopeTime t;
  this->doLoadData(m_trainFile, m_TrainData);
  this->doLoadData(m_testFile, m_PredictData);
  this->doLoadAnswer();

  std::cerr << "* TrainData: (" << m_TrainData.size() << ", "
            << m_TrainData[0].features.size() << ")\n";
  std::cerr << "* TestData: (" << m_PredictData.size() << ", "
            << m_PredictData[0].features.size() << ")\n";
  std::cerr << "* 正样本: " << m_positiveSamples << "\n";
  std::cerr << "* 负样本: " << m_negativeSamples << "\n";
  t.LogTime();
  std::cerr << "--------------------------------------\n";
}

void Simulation::SaveAnswer() {
  std::ofstream fout(m_predictFile);
  assert(fout);
  for (auto &it : m_PredictData) {
    fout << it.label << "\n";
  }
  fout.close();
}

int main() {
  std::cout << std::fixed << std::setprecision(4);
  std::cerr << std::fixed << std::setprecision(4);

  ScopeTime t;
  Simulation *simulation = new Simulation();

  simulation->LoadData();
  simulation->Train();
  simulation->Predict();
  simulation->SaveAnswer();
  simulation->Score();

  t.LogTime();
  delete simulation;
  return 0;
}