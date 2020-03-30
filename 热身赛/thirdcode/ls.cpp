#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#define NUM_THREADS 8

using namespace std;

//矩阵操作
struct Matrix {
  typedef vector<float> Mat1;
  typedef vector<vector<float>> Mat2;
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
static float Dot(const Matrix::Mat1 &mat1, const Matrix::Mat1 &mat2) {
  int n = mat1.size();
  float ans = 0;
  for (int i = 0; i < n; i += 16) {
    for (int j = 0; j < 16; j++) {
      ans += mat1[i + j] * mat2[i + j];
    }
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
  Matrix::Mat1 mat(n);
  for (int i = 0; i < n; i++) {
    mat[i] = mat1[i] - mat2[i];
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
  inline int predict(Matrix::Mat1 &vt);  //预测
  void judge();                          //判分
  void getNumber(int start);             //获取训练、测试数据信息

 public:
  // vector<char> answer;  //进行比对的答案
  char answer[40000];
  string trainFile;
  string testFile;
  string answerFile;
  string predictOutFile;

 public:
  int testLineSize = 6000;
  int feature = 400;    //特征数
  int trainNum = 1580;  //参与训练样本数
  float totTime = 0;    //运行时间
  int predictNum = 100000;
  int featureId;
  int testNum = 6000;
  float rho = 0;
  float pLabel0;
  float pLabel1;
  int items;
  vector<vector<float>> mu[2];
  vector<vector<int>> trainLabel;
  vector<vector<float>> delta;
  vector<vector<pair<int, int>>> threadStart;
  float pw[10][5] = {0, 0,   0,    0,     0,      1, 0.1, 0.01, 0.001, 0.0001,
                     2, 0.2, 0.02, 0.002, 0.0002, 3, 0.3, 0.03, 0.003, 0.0003,
                     4, 0.4, 0.04, 0.004, 0.0004, 5, 0.5, 0.05, 0.005, 0.0005,
                     6, 0.6, 0.06, 0.006, 0.0006, 7, 0.7, 0.07, 0.007, 0.0007,
                     8, 0.8, 0.08, 0.008, 0.0008, 9, 0.9, 0.09, 0.009, 0.0009};
  float pwChar[255][5];

 public:
  void writeData(Matrix::Mat1 &vt, string &file);
  void loadTrainData();
  void loadTestData(const string &file, int &lineSize, int pid);
  void LoadChar(const vector<char> &ch);
  void train();
  void normalization(int &pid, int type);
  void getReadId(char *buffer);
  void threadLoadData(char *buffer, int pid);
  void threadPredict(char *buffer, int pid, int start, int end, int lineSize);
  void LoadChar(char *buffer, int &pid, int &start, int &end);
};

LR::LR(string trainF, string testF, string predictOutF, string answerF) {
  trainFile = trainF;
  testFile = testF;
  predictOutFile = predictOutF;
  answerFile = answerF;
  init();
}

void LR::getReadId(char *buffer) {
  int now = 0, pre, threadId = 0, j = 0,
      p = (trainNum + NUM_THREADS - 1) / NUM_THREADS, circle = 0;

  for (int i = 0; i < NUM_THREADS; i++) {
    if (now + p <= trainNum) {
      threadStart.emplace_back(vector<pair<int, int>>(p));
      now += p;
    } else {
      threadStart.emplace_back(vector<pair<int, int>>(trainNum - now));
    }
  }
  now = 0;
  for (int i = 0; i < trainNum; ++i) {
    pre = now;
    now += testLineSize;
    while (buffer[now] != '\n') ++now;
    threadStart[threadId][j++] = make_pair(pre, now - 1);
    circle++;
    if (circle == p) {
      threadId++;
      j = 0;
      circle = 0;
    }
    now++;
  }
}
void LR::init() {
  // cout << "开始初始化\n";
  // cout << "NUM_THREADS: " << NUM_THREADS << "\n";
  // int start = clock();
  items = 64 / sizeof(float);
  for (int i = 0; i < NUM_THREADS; i++) {
    mu[0].emplace_back(vector<float>(feature, 0));
    mu[1].emplace_back(vector<float>(feature, 0));
    delta.emplace_back(vector<float>(2, 0));
    trainLabel.emplace_back(vector<int>(2, 0));
  }

  for (int i = '0'; i <= '9'; i++) {
    for (int j = 0; j < 5; j++) {
      pwChar[i][j] = pw[i - '0'][j];
    }
  }
  featureId = feature * 6;

  loadTrainData();

  // int mid = clock();
  // cout << "Load Train Data time: " << (float)(mid - start) / CLOCKS_PER_SEC
  //      << " s\n";

  train();
  // int mid2 = clock();
  // cout << "Train Time: " << (float)(mid2 - mid) / CLOCKS_PER_SEC << " s\n";

  // loadTestData(testFile, testLineSize);

  // mid = clock();
  // cout << "Predict Time: " << (float)(mid - mid2) / CLOCKS_PER_SEC <<
  // "s\n"; int end = clock(); cout << " Running Time: " << (float)(end -
  // start) / CLOCKS_PER_SEC
  //      << " s\n";
}

//////////////////

void LR::LoadChar(char *buffer, int &pid, int &start, int &end) {
  int now = start, id = 0, r = 0;
  float num = 0, sum = 0;
  bool flag = false;

  int type = buffer[end] - '0';
  trainLabel[pid][type]++;
  while (id < feature) {
    if (buffer[now] == '-') {
      now++;
      num = pwChar[buffer[now]][0] + pwChar[buffer[now + 2]][1] +
            pwChar[buffer[now + 3]][2] + pwChar[buffer[now + 4]][3];
      mu[type][pid][id++] -= num;
      sum -= num;
    } else {
      num = pwChar[buffer[now]][0] + pwChar[buffer[now + 2]][1] +
            pwChar[buffer[now + 3]][2] + pwChar[buffer[now + 4]][3];
      mu[type][pid][id++] += num;
      sum += num;
    }
    now += 6;
  }
  delta[pid][type] += sum;
}

void LR::threadLoadData(char *buffer, int pid) {
  for (auto &x : threadStart[pid]) {
    LoadChar(buffer, pid, x.first, x.second);
  }
}

void LR::loadTrainData() {
  struct stat sb;
  int fd = open(trainFile.c_str(), O_RDONLY);
  fstat(fd, &sb);
  char *buffer = (char *)mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  getReadId(buffer);
  vector<thread> td(NUM_THREADS);
  for (int i = 0; i < NUM_THREADS; i++) {
    td[i] = thread(&LR::threadLoadData, this, buffer, i);
  }
  for (auto &t : td) {
    t.join();
  }
  close(fd);
}

void LR::threadPredict(char *buffer, int pid, int start, int end,
                       int lineSize) {
  int id, initId = start * lineSize, nowId, up, now = 0, j, k, r;
  float sum;
  cout << pid << " " << start << " " << end << "\n";
  for (int i = start; i < end; i += items) {
    up = items;
    if (i + items > end) up = end - i;
    for (k = 0; k < up; k++) {
      id = 0;
      sum = pLabel1;
      for (j = 0; j < featureId; j += 60) {
        nowId = initId + j;
        for (r = 0; r < 60; r += 6) {
          sum +=
              (pwChar[buffer[nowId + r]][0] + pwChar[buffer[nowId + r + 2]][1] +
               pwChar[buffer[nowId + r + 3]][2]) *
              mu[0][0][id];
          ++id;
        }
      }
      answer[(i + k) << 1 | 1] = '\n';
      answer[(i + k) << 1] = sum > 0 ? '1' : '0';
      initId += lineSize;
    }
  }
}

void LR::loadTestData(const string &file, int &lineSize, int pid) {
  struct stat sb;
  int fd = open(file.c_str(), O_RDONLY);
  fstat(fd, &sb);
  char *buffer = (char *)mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);
  int linenum = sb.st_size / lineSize;
  int pre = 0, line = (linenum + NUM_THREADS - 1) / NUM_THREADS;
  // answer = vector<char>(linenum * 2);
  pre = line * pid;
  threadPredict(buffer, pid, pre, min(pre + line, linenum), lineSize);
}

////////////////////
// void LR::normalization(int &pid, int type) {
//     float sum = 0;
//     trainLabel[pid][type]++;
//     for (int i = 0; i < features[pid].size(); i += items) {
//         for (int j = 0; j < items; j++) {
//             mu[type][pid][i + j] += features[pid][i + j];
//             sum += features[pid][i + j];
//         }
//     }
//     delta[pid][type] += sum;
// }

void LR::train() {
  for (int i = 1; i < NUM_THREADS; i++) {
    trainLabel[0][0] += trainLabel[i][0];
    delta[0][0] += delta[i][0];
    delta[0][1] += delta[i][1];
  }
  delta[0][0] += 1;
  delta[0][1] += 1;

  // printf("delta: %.6f %.6f\n", delta[0][0], delta[0][1]);

  pLabel0 = 1.0 * trainLabel[0][0] / trainNum;
  pLabel1 = 1.0 - pLabel0;
  pLabel0 = log(pLabel0);
  pLabel1 = log(pLabel1);

  // cout << "log p: " << pLabel0 << " " << pLabel1 << "\n";

  float al0 = log(1.0 / delta[0][0]), al1 = log(1.0 / delta[0][1]);
  int nowId;
  for (int i = 0; i < feature; i += items) {
    for (int k = 0; k < items; k++) {
      nowId = i + k;
      for (int j = 1; j < NUM_THREADS; j++) {
        mu[0][0][nowId] += mu[0][j][nowId];
        mu[1][0][nowId] += mu[1][j][nowId];
      }
      // mu[0][0][nowId] += 1;
      // mu[1][0][nowId] += 1;

      mu[0][0][nowId] = log(mu[0][0][nowId] + 1) + al0;
      mu[1][0][nowId] = log(mu[1][0][nowId] + 1) + al1;
    }
  }
  mu[0][0] = mu[1][0] - mu[0][0];
  pLabel1 -= pLabel0;
}

void LR::judge() {
  vector<char> answer2;
  char x;
  int cor = 0, id = 0;
  ifstream fin(answerFile);
  while (fin) {
    fin >> x;
    answer2.emplace_back(x);
  }
  fin.close();
  cout << answer2.size() << "\n";
  cout << answer[20000] << "\n";
  for (int i = 0; i < 10000; i++) {
    if (answer[i * 2] == answer2[i]) cor++;
    // else {
    //     cout << answer[i] << " " << answer2[i] << "\n";
    // }
  }
  cout << "准确率: " << 1.0 * cor / 40000 << "\n";
}
int main(int argc, char *argv[]) {
  // srand((unsigned)time(NULL));

  string trainFile = "../data/train_data.txt";
  string testFile = "../data/test_data.txt";
  string predictFile = "../data/result.txt";
  string answerFile = "../data/answer.txt";

  // string trainFile = "/data/train_data.txt";
  // string testFile = "/data/test_data.txt";
  // string predictFile = "/projects/student/result.txt";
  // string answerFile = "/projects/student/answer.txt";

  LR logist(trainFile, testFile, predictFile, answerFile);

  pid_t child1 = 0, child2 = 0, child3 = 0;
  pid_t child4 = 0, child5 = 0, child6 = 0, child7 = 0;

  child1 = fork();
  if (child1) child2 = fork();
  if (child2) child3 = fork();
  if (child3) child4 = fork();
  if (child4) child5 = fork();
  if (child5) child6 = fork();
  if (child6) child7 = fork();

  int pid = 0;
  if (!child1) {
    pid = 1;
  } else if (!child2) {
    pid = 2;
  } else if (!child3) {
    pid = 3;
  } else if (!child4) {
    pid = 4;
  } else if (!child5) {
    pid = 5;
  } else if (!child6) {
    pid = 6;
  } else if (!child7) {
    pid = 7;
  } else {
    pid = 0;
  }

  if (pid <= 7) {
    logist.loadTestData(logist.testFile, logist.testLineSize, pid);
    if (pid) exit(0);
  }

  waitpid(child1, NULL, 0);
  waitpid(child2, NULL, 0);
  waitpid(child3, NULL, 0);
  waitpid(child4, NULL, 0);
  waitpid(child5, NULL, 0);
  waitpid(child6, NULL, 0);
  waitpid(child7, NULL, 0);

  for (int i = 0; i < 40000; i++) {
    cout << i << " " << logist.answer[i] << "\n";
  }
  // FILE *fp = fopen(logist.predictOutFile.c_str(), "w");
  // fwrite(logist.answer, 40000, 1, fp);
  // fclose(fp);

  logist.judge();
  return 0;
}