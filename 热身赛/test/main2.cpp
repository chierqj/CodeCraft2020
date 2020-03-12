// #include <Windows.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <thread>
#include <vector>

// #include "Ctime.h"
// CTimer timer;

#define TEST
#define thnum 8
#define linenumi 20000
using namespace std;

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

struct Data {
  vector<double> features;
  int label;
  Data(vector<double> f, int l) : features(f), label(l) {}
};

struct Param {};

class testway {
 public:
  void train();
  void predict();
  int loadModel();
  int storeModel();
  testway(string trainFile, string testFile, string predictOutFile);

 private:
  vector<Data> trainDataSet[thnum];
  vector<Data> testDataSet;
  vector<int> predictVec;
  Param param;
  string trainFile;
  string testFile;
  string predictOutFile;
  string weightParamFile = "modelweight.txt";

 private:
  bool init();
  void initParam();
  int calmn();
  bool multiload(int threadnum, int linenum, int linesizei);
  bool TrainDatajoin();

 private:
  int featuresNum;
  int linesize;
  int linenum;
  const double miuInitV = 0;
  const double sigmaInitV = 0;
  const double stepSize = 0.1;
  const int maxIterTimes = 3;
  const double predictTrueThresh = 1.;
  const int train_show_step = 10;
};

testway::testway(string trainF, string testF, string predictOutF) {
  trainFile = trainF;
  testFile = testF;
  predictOutFile = predictOutF;
  featuresNum = 0;

  calmn();  //ͳ���ļ�����������
  init();
}

void testway::initParam() {}

bool testway::init() {
  for (int i = 0; i < thnum; i++) {
    trainDataSet[thnum].clear();
  }

  bool status = TrainDatajoin();
  //	loadTrainData();
  if (status != true) {
    return false;
  }
  featuresNum = trainDataSet[0][0].features.size();
  ////param.NormalSet.clear();

  initParam();
  return true;
}

bool testway::multiload(int threadnum, int linenum, int linesizei) {
  ///////////////////////
  // int threadnum = 0;
  // int linesizei = 60002;
  // int linenumi = 1000;
  //����ƫ����
  int datamove = threadnum * linesizei * linenum;
  //////////////////////////
  char a[10];
  int c, i = 0;
  int n = 0, m = 0, mi = 0;
  vector<double> feature;
  feature.clear();

  FILE *pFile;
  int len;
  char buffer[500000];
  int bufsize = 500000;
  pFile = fopen(trainFile.c_str(), "r");
  fseek(pFile, datamove, SEEK_SET);
  if (pFile == NULL) {
    return 0;
  } else {
    do {
      len = fread(buffer, 1, bufsize, pFile);
      for (int j = 0; j < len; j++) {
        mi++;
        if (buffer[j] == '\n') {
          n++;
          if (m == 0) m = mi;
          i = 0;
          int temp = a[0] - '0';
          trainDataSet[threadnum].push_back(Data(feature, temp));
          if (n >= linenum) {
            return 1;
          }
          feature.clear();
        } else {
          if (buffer[j] == ',') {
            a[i] = '\n';
            double temp = strtod(a, NULL);
            feature.push_back(temp);
            i = 0;
          } else {
            a[i] = buffer[j];
            i++;
          }
        }
      }

    } while (!feof(pFile));
  }

  fclose(pFile);
  return 1;
}
bool testway::TrainDatajoin() {
  vector<std::thread> vec_threads;
  for (int i = 0; i < thnum; ++i) {
    int linenumtemp = linenum / thnum;

    std::thread th(&testway::multiload, this, i, linenumi, linesize);
    vec_threads.emplace_back(std::move(th));  // push_back() is also OK
  }

  auto it = vec_threads.begin();
  for (; it != vec_threads.end(); ++it) {
    (*it).join();
  }
  for (int i = 1; i < thnum; ++i) {
    trainDataSet[0].insert(trainDataSet[0].end(), trainDataSet[i].begin() + 1,
                           trainDataSet[i].end() - 1);
  }
  return 1;
}
int testway::calmn() {
  int i = 0;
  int n = 0, m = 0, mi = 0;

  FILE *pFile;
  int len;
  char buffer[500000];
  int bufsize = 500000;
  pFile = fopen(trainFile.c_str(), "r");
  if (pFile == NULL) {
    return 0;
  } else {
    do {
      len = fread(buffer, 1, bufsize, pFile);
      for (int j = 0; j < len; j++) {
        mi++;
        if (buffer[j] == '\n') {
          n++;
          if (m == 0) {
            // fclose(pFile);
            m = mi;
            linesize = mi;
            // return 1;
          }
        }
      }
    } while (!feof(pFile));
  }

  linesize = m;
  linenum = n;
  cout << "(" << n << ", " << m << ")" << endl;
  fclose(pFile);
  // cout << "�����С��ʱ��" << timer.time_out() << endl;
  // timer.time_in();
  return 1;
}
int main(int argc, char *argv[]) {
  // timer.time_in();
  ScopeTime t;
  vector<int> answerVec;
  vector<int> predictVec;
  int correctCount;
  double accurate;
  string trainFile = "../code/data/train_data.txt";
  string testFile = "../code/data/test_data.txt";
  string predictFile = "../code/data/result.txt";

  string answerFile = "../code/data/answer.txt";

  testway testcal(trainFile, testFile, predictFile);

  // cout << "���ݶ�ȡʱ�䣺" << timer.time_out() << endl;
  t.LogTime();
  // timer.time_in();

  getchar();
  return 0;
}
