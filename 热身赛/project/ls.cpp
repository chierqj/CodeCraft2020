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
#include <cstdlib>
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

using namespace std;

float pw[10];
//矩阵操作
struct Matrix {
  typedef vector<float> Mat1;
  typedef vector<vector<float> > Mat2;
};
void out(Matrix::Mat1 mat) {
  std::cout << "{";
  for (auto &x : mat) {
    cout << x << " ";
  }
  cout << "}\n";
}
void out(Matrix::Mat2 mat) {
  std::cout << "[";
  for (auto &x : mat) {
    out(x);
  }
  std::cout << "]\n";
}
//点乘
static float Dot(const Matrix::Mat1 &mat1, const Matrix::Mat1 &mat2) {
  int n = mat1.size();
  float ans = 0;
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
  int predict(Matrix::Mat1 &vt);  //预测
  void judge();                   //判分
  void getNumber(int start);      //获取训练、测试数据信息

 private:
  Matrix::Mat2 trainX;  //训练数据集
  Matrix::Mat1 trainY;  //训练集标签
  Matrix::Mat2 testX;   //测试数据集
  Matrix::Mat1 testY;   //预测的结果
  vector<int> answer;   //进行比对的答案
  Matrix::Mat1 pred;    //预测的sigmod结果
  string trainFile;
  string testFile;
  string answerFile;
  string predictOutFile;

 private:
  int testLineSize = 6000;
  int feature = 100;    //特征数
  int trainNum = 3000;  //参与训练样本数
  int YLabTrain = 0;    //测试结果中正标签
  int NLabTrain = 0;    //测试结果中负标签
  int YLabTest = 0;     //测试结果中正标签
  int NLabTest = 0;     //测试结果中负标签
  float totTime = 0;    //运行时间
  float rho = 0;
  float pLabel0;
  float pLabel1;
  vector<float> mu[2];
  vector<int> trainLabel;
  vector<float> delta;

 private:
  void writeData(Matrix::Mat1 &vt, string &file);
  void loadTrainData(const string &file);
  void loadTestData(const string &file, int &lineSize);
  void LoadChar(const vector<char> &ch);
  void train();
  void normalization(Matrix::Mat1 &vt, int type);
};

inline void ShuffleVector(vector<int> &vt) {
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  shuffle(vt.begin(), vt.end(), std::default_random_engine(seed));
}

LR::LR(string trainF, string testF, string predictOutF, string answerF) {
  trainFile = trainF;
  testFile = testF;
  predictOutFile = predictOutF;
  answerFile = answerF;
  init();
}

void LR::init() {
  cout << "开始初始化\n";
  int start = clock();
  mu[0] = vector<float>(feature, 1);
  mu[1] = vector<float>(feature, 1);
  delta = vector<float>(2, 1);
  trainLabel = vector<int>(2, 0);
  loadTrainData(trainFile);
  train();
  loadTestData(testFile, testLineSize);
  int end = clock();
  cout << "Running Time: " << (float)(end - start) / CLOCKS_PER_SEC << " s\n";
}

//////////////////

// double Todouble(string &s) {
//     int n = s.size(), i = n - 1;
//     double ans = 0;
//     if (s[0] == '-') {
//         return 0;
//     } else if (s[0] == '0') {
//         while (i > 1) {
//             ans = ans + (s[i--] - '0');
//             ans /= 10;
//         }
//         return ans;
//     } else {
//         return 1;
//     }
// }

// void LR::loadTrainData(const string &file) {
//     ifstream fin(file);
//     vector<float> features;
//     string line, temp;
//     int label, num = 0;
//     while (fin) {
//         getline(fin, line);
//         if (num == trainNum) {
//             break;
//         }
//         if (line.empty()) {
//             continue;
//         }
//         istringstream iss(line);
//         features.clear();
//         while (getline(iss, temp, ',')) {
//             features.emplace_back(Todouble(temp));
//         }
//         label = features.back();
//         features.pop_back();

//         normalization(features, label);
//         features.clear();
//         num++;
//     }
//     fin.close();
// }

void LR::loadTrainData(const string &file) {
  struct stat sb;
  int fd = open(file.c_str(), O_RDONLY);
  fstat(fd, &sb);
  char *buffer = (char *)mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  int r, id, now = 0, p = 0, label;
  bool flag;
  float num, sum;
  Matrix::Mat1 features(feature);
  ofstream fout(predictOutFile);
  for (int i = 0; i < trainNum; ++i) {
    id = 0;
    num = 0;
    r = 0;
    flag = false;
    sum = 0;
    while (id < feature) {
      if (buffer[now] == ',') {
        if (flag) {
          features[id] = 0 - num;
        } else {
          features[id] = num;
        }
        num = 0;
        id++;
        r = 0;
        flag = false;
      } else if (buffer[now] == '-') {
        flag = true;
      } else if (buffer[now] != '.') {
        p = (buffer[now] - '0');
        num += pw[r] * p;
        r++;
      }
      now++;
    }
    while (buffer[now + 1] != '\n') ++now;
    label = buffer[now] - '0';
    now += 2;
    normalization(features, label);
  }
}

void LR::loadTestData(const string &file, int &lineSize) {
  struct stat sb;
  int fd = open(file.c_str(), O_RDONLY);
  fstat(fd, &sb);
  char *buffer = (char *)mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  long long linenum = sb.st_size / lineSize;
  int x1, x2, x3, x4, id;
  long long initId;
  float num;
  Matrix::Mat1 features(feature);
  ofstream fout(predictOutFile);
  for (int i = 0; i < linenum; ++i) {
    id = 0;
    initId = (long long)i * lineSize;
    for (int j = 0; j < lineSize; j += 6) {
      x1 = buffer[initId + j] - '0';
      x2 = buffer[initId + j + 2] - '0';
      x3 = buffer[initId + j + 3] - '0';
      x4 = buffer[initId + j + 4] - '0';
      if (id < feature)
        // features[id++] = x1 + (float)(x2 * 100 + x3 * 10 + x4) /
        // 1000;
        features[id++] = x1 * 1.0 + x2 * 0.1 + x3 * 0.01 + x4 * 0.001;
      else
        break;
    }
    fout << predict(features) << "\n";
  }
  fout.close();
}

////////////////////
void LR::normalization(Matrix::Mat1 &vt, int type) {
  int id = 0;
  trainLabel[type]++;
  float sum = 0;
  for (auto &x : vt) {
    mu[type][id] += x;
    sum += x;
    id++;
  }
  delta[type] += sum;
}

void LR::train() {
  pLabel0 = 1.0 * trainLabel[0] / trainNum;
  pLabel1 = 1.0 - pLabel0;
  pLabel0 = log(pLabel0);
  pLabel1 = log(pLabel1);
  float al0 = 1.0 / delta[0], al1 = 1.0 / delta[1];
  for (int i = 0; i < feature; i++) {
    mu[0][i] = log(mu[0][i] * al0);
    mu[1][i] = log(mu[1][i] * al1);
  }
}
int LR::predict(Matrix::Mat1 &vt) {
  float sum0, sum1;
  sum0 = Dot(vt, mu[0]);
  sum1 = Dot(vt, mu[1]);
  sum0 += pLabel0;
  sum1 += pLabel1;
  // std::cout << sum0 << ", " << sum1 << "\n";
  if (sum1 > sum0) return 1;
  return 0;
}

void LR::judge() {
  vector<int> answer, result;
  int x, cor = 0;
  ifstream fin(answerFile);
  while (fin) {
    fin >> x;
    answer.emplace_back(x);
  }
  fin.close();
  ifstream fin2(predictOutFile);
  while (fin2) {
    fin2 >> x;
    result.emplace_back(x);
  }
  fin2.close();
  for (int i = 0; i < answer.size(); i++) {
    if (answer[i] == result[i]) cor++;
  }
  cout << "准确率: " << 1.0 * cor / answer.size() << "\n";
}
int main(int argc, char *argv[]) {
  std::cout << std::fixed << std::setprecision(3);

  // srand((unsigned)time(NULL));
  pw[0] = 1;
  for (int i = 1; i < 8; i++) {
    pw[i] = pw[i - 1] * 0.1;
  }
  string trainFile = "../data/train_data.txt";
  string testFile = "../data/test_data.txt";
  string predictFile = "../data/result.txt";
  string answerFile = "../data/answer.txt";

  // string trainFile = "/data/train_data.txt";
  // string testFile = "/data/test_data.txt";
  // string predictFile = "/projects/student/result.txt";
  // string answerFile = "/projects/student/answer.txt";

  LR logist(trainFile, testFile, predictFile, answerFile);
  // logist.getNumber(start_time);
  logist.judge();
  return 0;
}