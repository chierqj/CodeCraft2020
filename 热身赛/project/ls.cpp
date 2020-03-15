#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <thread>
#include <vector>
using namespace std;

int pw[10];
//矩阵操作
struct Matrix {
  typedef vector<double> Mat1;
  typedef vector<vector<double> > Mat2;
};
void out(Matrix::Mat1 mat) {
  for (auto &x : mat) {
    cout << x << " ";
  }
  cout << "\n";
}
void out(Matrix::Mat2 mat) {
  for (auto &x : mat) {
    out(x);
  }
}
//点乘
static double Dot(const Matrix::Mat1 &mat1, const Matrix::Mat1 &mat2) {
  int n = mat1.size();
  double ans = 0;
  for (int i = 0; i < n; i++) {
    ans += mat1[i] * mat2[i];
  }
  return ans;
};
//乘法
static Matrix::Mat1 operator*(const Matrix::Mat2 &mat1,
                              const Matrix::Mat1 &mat2) {
  int n = mat1.size();
  Matrix::Mat1 mat;
  for (const auto &x : mat1) {
    mat.emplace_back(Dot(x, mat2));
  }
  return mat;
};
static Matrix::Mat2 operator*(const Matrix::Mat1 &mat1,
                              const Matrix::Mat1 &mat2) {
  int n = mat1.size(), m = mat2.size(), id = 0;
  Matrix::Mat2 mat(n);
  for (const auto &x : mat1) {
    for (const auto &y : mat2) {
      mat[id].emplace_back(x * y);
    }
    id++;
  }
  return mat;
};
//转置
static Matrix::Mat2 T(const Matrix::Mat2 &mat1) {
  int n = mat1.size(), m = mat1[0].size();
  Matrix::Mat2 mat(m, Matrix::Mat1(n));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      mat[j][i] = mat1[i][j];
    }
  }
  return mat;
};
static Matrix::Mat1 operator-(const Matrix::Mat1 &mat1,
                              const Matrix::Mat1 &mat2) {
  int n = mat1.size();
  Matrix::Mat1 mat;
  for (int i = 0; i < n; i++) {
    mat.emplace_back(mat1[i] - mat2[i]);
  }
  return mat;
}
static Matrix::Mat2 operator*(const Matrix::Mat2 &mat1,
                              const Matrix::Mat2 &mat2) {
  Matrix::Mat2 mat2T = T(mat2);
  int n = mat1.size(), m = mat2[0].size();
  Matrix::Mat2 mat(n, Matrix::Mat1(m));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      mat[i][j] = Dot(mat1[i], mat2T[j]);
    }
  }
  return mat;
};

class LR {
 public:
  LR(string trainF, string testF, string predictOutF, string answerF);
  void init();
  void predict();             //预测
  void judge(string &file);   //判分
  void getNumber(int start);  //获取训练、测试数据信息

 private:
  Matrix::Mat2 trainX;  //训练数据集
  Matrix::Mat1 trainY;  //训练集标签
  Matrix::Mat2 testX;   //测试数据集
  Matrix::Mat1 testY;   //预测的结果
  vector<int> answer;   //进行比对的答案
  Matrix::Mat1 pred;    //预测的sigmod结果
  Matrix::Mat1 weight;  //权重矩阵
  string trainFile;
  string testFile;
  string answerFile;
  string predictOutFile;

 private:
  double learningRate = 0.01;  //学习率
  int iterations = 130;        //训练轮数
  int feature = 0;             //特征数
  int trainNum = 800;          //样本数
  int YLabTrain = 0;           //样本中正标签
  int NLabTrain = 0;           //样本中负标签
  int YLabTest = 0;            //测试结果中正标签
  int NLabTest = 0;            //测试结果中负标签
  double totTime = 0;          //运行时间
  char *m_Buffer;
  std::vector<Matrix::Mat2> m_ThreadData;

 private:
  void loadData(const string &file, Matrix::Mat2 &vt1, Matrix::Mat1 &vt2);
  void writeWeight(Matrix::Mat1 &vt, string &file);
  void writeData(Matrix::Mat1 &vt, string &file);
  void LoadChar(const vector<char> &c, Matrix::Mat2 &vt1, Matrix::Mat1 &vt2);
  void changeWeight(Matrix::Mat1 &vt);
  void NewloadData(const string &file, Matrix::Mat2 &vt1, Matrix::Mat1 &vt2);
  void init_weight();
  void train();
  int getLab(double &z);
  double sigmod(double &x);
  void handleThread(int pid, unsigned int left, unsigned int right);
};

double LR::sigmod(double &x) { return 1.0 / (1.0 + exp(-x)); }
LR::LR(string trainF, string testF, string predictOutF, string answerF) {
  trainFile = trainF;
  testFile = testF;
  predictOutFile = predictOutF;
  answerFile = answerF;
  init();
  train();
}

void LR::getNumber(int start) {
  YLabTest = 0;
  YLabTrain = 0;
  NLabTest = 0;
  NLabTrain = 0;
  cout << "================\n";
  for (auto &x : trainY) {
    if (x == 1) {
      YLabTrain++;
    } else
      NLabTrain++;
  }
  cout << "================\n";
  for (auto &x : testY) {
    if (x > 0.5) {
      YLabTest++;
    } else
      NLabTest++;
  }
  int end = clock();
  totTime = (double)(end - start) / CLOCKS_PER_SEC;
  cout << "================\n";
  cout << "训练集样本数：" << trainNum << "\n";
  cout << "训练集特征数：" << feature << "\n";
  cout << "训练集正样本数：" << YLabTrain << "\n";
  cout << "训练集负样本数：" << NLabTrain << "\n";
  cout << "================\n";
  cout << "测试集样本数：" << testY.size() << "\n";
  cout << "测试集正样本数：" << YLabTest << "\n";
  cout << "测试集负样本数：" << NLabTest << "\n";
  cout << "运行时间：" << totTime << "s\n";
  cout << "================\n";
}

void LR::init_weight() {
  weight = Matrix::Mat1(feature, 1);
  weight.emplace_back(0);
  return;
  int limit = sqrt(1.0 / feature) * 10000;
  double p;
  for (int i = 0; i < feature; i++) {
    p = (rand() % (limit * 2) - limit) * 1.0 / 10000;
    weight.emplace_back(p);
  }
  weight.emplace_back(0);
}

void LR::init() {
  // cout << "开始初始化\n";
  NewloadData(trainFile, trainX, trainY);
  loadData(testFile, testX, testY);
  // cout << trainX.size() << " " << testX.size() << "\n";
  feature = trainX[0].size() - 1;
  init_weight();
}

//////////////////

double Todouble(string &s) {
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
void LR::NewloadData(const string &file, Matrix::Mat2 &vt1, Matrix::Mat1 &vt2) {
  ifstream fin(file);
  vector<double> features;
  string line, temp;
  int label;
  while (fin) {
    getline(fin, line);
    if (vt1.size() == trainNum) {
      break;
    }
    if (line.empty()) {
      continue;
    }
    istringstream iss(line);
    features.clear();
    while (getline(iss, temp, ',')) {
      features.emplace_back(Todouble(temp));
    }
    label = features.back();
    features.pop_back();
    features.emplace_back(1);
    vt1.emplace_back(features);
    vt2.emplace_back(label);
    features.clear();
  }
  fin.close();
}

////////////////////

void LR::LoadChar(const vector<char> &ch, Matrix::Mat2 &vt1,
                  Matrix::Mat1 &vt2) {
  bool ne = false, dot = false, one = false;
  int r = 1;
  double sum = 0, x;
  Matrix::Mat1 features;
  for (auto &c : ch) {
    if (c == ',' || c == '\n') {
      if (ne) {
        sum = 0;
        ne = false;
      }
      if (one) {
        sum = 1;
        one = false;
      }
      dot = false;
      r = 1;
      if (c == ',') {
        features.emplace_back(sum);
      } else {
        features.emplace_back(sum);
        features.emplace_back(1);
        vt1.emplace_back(features);
        features.clear();
      }
      sum = 0;
    } else {
      if (ne || one) {
        continue;
      }
      if (c == '-') {
        ne = true;
        continue;
      } else if (c == '.') {
        dot = true;
        continue;
      }
      x = c - '0';
      if (!dot) {
        if (c == '0')
          sum = 0;
        else
          one = true;

      } else {
        sum += x / pw[r];
        r++;
      }
    }
  }
}
void LR::handleThread(int pid, unsigned int left, unsigned int right) {
  Matrix::Mat1 features;
  for (unsigned int i = left; i < right; i += 6) {
    double x1 = m_Buffer[i] - '0';
    double x2 = m_Buffer[i + 2] - '0';
    double x3 = m_Buffer[i + 3] - '0';
    double x4 = m_Buffer[i + 4] - '0';
    double num = x1 + x2 / 10 + x3 / 100 + x4 / 1000;
    features.emplace_back(num);
    if (m_Buffer[i + 5] == '\n') {
      m_ThreadData[pid].emplace_back(features);
      features.clear();
    }
  }
}
void LR::loadData(const string &file, Matrix::Mat2 &vt1, Matrix::Mat1 &vt2) {
  ifstream fin(file, ios::binary);
  unsigned int bufsize = fin.seekg(0, ios::end).tellg();
  m_Buffer = new char[bufsize];
  fin.seekg(0, ios::beg).read(m_Buffer, static_cast<streamsize>(bufsize));
  fin.close();

  Matrix::Mat1 features;
  int p = 0;
  int linesize = 0;
  while (p < bufsize) {
    double x1 = m_Buffer[p] - '0';
    double x2 = m_Buffer[p + 2] - '0';
    double x3 = m_Buffer[p + 3] - '0';
    double x4 = m_Buffer[p + 4] - '0';
    double num = x1 + x2 / 10 + x3 / 100 + x4 / 1000;
    features.emplace_back(num);
    if (m_Buffer[p + 5] == '\n') {
      vt1.emplace_back(features);
      features.clear();
      linesize = p + 6;
      break;
    }
    p += 6;
  }

  int nthread = 4;
  m_ThreadData.resize(nthread);

  int linenum = bufsize / linesize;
  int blocksize = linenum / nthread;
  int start = 1;
  std::vector<std::thread> Thread(nthread);
  for (int i = 0; i < nthread; i++) {
    int end;
    if (i == nthread - 1) {
      end = linenum;
    } else {
      end = start + blocksize;
    }
    unsigned int l = start * linesize, r = end * linesize;
    Thread[i] = std::thread(&LR::handleThread, this, i, l, r);
    start += blocksize;
  }
  for (auto &it : Thread) {
    it.join();
  }
  for (auto &it : m_ThreadData) {
    vt1.insert(vt1.end(), it.begin(), it.end());
  }
}

void LR::writeWeight(Matrix::Mat1 &vt, string &file) {
  string line;
  ofstream fout(file);
  for (auto &x : vt) {
    fout << x << "\n";
  }
  fout.close();
}

void LR::train() {
  // cout << "================\n";
  // cout << "------开始训练------\n";
  Matrix::Mat2 mat2 = T(trainX);
  Matrix::Mat1 mat;
  Matrix::Mat1 sig;
  for (int i = 0; i < iterations; i++) {
    // if (i % 20 == 0) cout << "第 " << i << " 轮迭代\n";
    learningRate = 5.0 / (i + 1) + 0.01;
    mat = trainX * weight;
    for (auto &x : mat) {
      x = sigmod(x);
    }
    sig = mat - trainY;
    sig = mat2 * sig;
    for (int j = 0; j <= feature; j++) {
      weight[j] -= sig[j] * learningRate;
    }
  }
  // string modelweight = "modelweight.txt";
  // writeWeight(weight, modelweight);
  // cout << "----------------\n";
  // cout << "----------------\n";
  // cout << "------训练结束------\n";
  // cout << "================\n";
}

void LR::writeData(Matrix::Mat1 &vt, string &file) {
  string line;
  ofstream fout(file);
  for (auto &x : vt) {
    fout << (x > 0.5 ? 1 : 0) << "\n";
  }
  fout.close();
}
void LR::judge(string &file) {
  ifstream fin(file);
  int x, cor = 0, n = testY.size();
  while (fin) {
    fin >> x;
    answer.emplace_back(x);
  }
  fin.close();
  for (int i = 0; i < testY.size(); i++) {
    if ((testY[i] > 0.5 ? 1 : 0) == answer[i]) cor++;
  }
  cout << "================\n";
  cout << "预测准确率为: " << cor * 1.0 / n << "\n";
  cout << "================\n";
}
void LR::predict() {
  testY = testX * weight;
  for (auto &x : testY) {
    x = sigmod(x);
  }
  writeData(testY, predictOutFile);
}

int main(int argc, char *argv[]) {
  srand((unsigned)time(NULL));
  pw[0] = 1;
  for (int i = 1; i < 5; i++) {
    pw[i] = pw[i - 1] * 10;
  }
  int start_time = clock();
  // string trainFile = "../data/train_data.txt";
  // string testFile = "../data/test_data.txt";
  // string predictFile = "../data/result.txt";
  // string answerFile = "../data/answer.txt";

  string trainFile = "/data/train_data.txt";
  string testFile = "/data/test_data.txt";
  string predictFile = "/projects/student/result.txt";
  string answerFile = "/projects/student/answer.txt";

  LR logist(trainFile, testFile, predictFile, answerFile);
  // cout << "================\n";
  // cout << "开始预测\n";
  logist.predict();

  // cout << "预测完成\n";
  // cout << "================\n";
  // logist.getNumber(start_time);
  // logist.judge(answerFile);
  return 0;
}