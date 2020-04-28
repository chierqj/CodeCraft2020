/*******
 * 多线程读 getline 读部分特征存入train_x再合并
 *
 * ********/
// #include <omp.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

using namespace std::chrono;
using namespace std;
#define THREADS_TRAIN 4
#define THREADS_TEST 4
#define FILE_SIZE 48217688
#define FEATRURENS_NUM 260
#define CacheLineSize 4

#define TEST1

#ifdef TEST
static const char *trainFile = "../data/train_data.txt";
static const char *predictFile = "../data/result.txt";
static const char *testFile = "../data/test_data.txt";
static const char *answerFile = "../data/answer.txt";

#else
static const char *trainFile = "/data/train_data.txt";
static const char *testFile = "/data/test_data.txt";
static const char *predictFile = "/projects/student/result.txt";
static const char *answerFile = "/projects/student/answer.txt";
#endif

int mean_fearture_1[THREADS_TRAIN][FEATRURENS_NUM] = {0};
int mean_fearture_0[THREADS_TRAIN][FEATRURENS_NUM] = {0};

int train_cnt[THREADS_TRAIN] = {0};
int label_1_arry[THREADS_TRAIN] = {0};

int test_x[THREADS_TEST][6000][FEATRURENS_NUM];
int test_cnt[THREADS_TEST] = {0};
#ifdef TEST
steady_clock::time_point start = steady_clock::now();
#endif
double get_accuary(const vector<int> &pred, const vector<int> &Y) {
  int correctCount = 0;
  size_t j = 0;
  for (; j < pred.size(); j += CacheLineSize) {
    for (; j < pred.size(); j++) {
      if (pred[j] == Y[j]) correctCount++;
    }
  }
  return ((double)correctCount) / pred.size();
}

int get_distance(const int *x1, const int *x0, const int val[]) {
  int sum = 0;
  size_t i = 0;
  for (i = 0; i < FEATRURENS_NUM; i += CacheLineSize) {
    sum += abs(x1[i] - val[i]) - abs(x0[i] - val[i]);
    sum += abs(x1[i + 1] - val[i + 1]) - abs(x0[i + 1] - val[i + 1]);
    sum += abs(x1[i + 2] - val[i + 2]) - abs(x0[i + 2] - val[i + 2]);
    sum += abs(x1[i + 3] - val[i + 3]) - abs(x0[i + 3] - val[i + 3]);
  }
  return sum < 300;
}

void thread_read_test(int id, long long start_point, long long volume) {
  ifstream is(testFile);
  is.seekg(start_point);
  string line;

  int line_count = 0;
  vector<vector<int>> thread_test;
  while (volume > 0 && getline(is, line)) {
    volume -= line.size() + 1;
    int count = 0;
    for (int j = 0; count < FEATRURENS_NUM;) {
      test_x[id][test_cnt[id]][count++] = (line[j + 2] - '0') * 100;
      j += 6;
    }
    test_cnt[id]++;
  }
  is.close();
#ifdef TEST
  cout << "test:" << id << " time: "
       << duration_cast<duration<double>>(steady_clock::now() - start).count()
       << endl;
#endif
}

void LoadTestLabels(vector<int> &test_y) {
  ifstream infile(answerFile);
  int i;
  while (infile >> i) test_y.push_back(i);
  infile.close();
}

void thread_read_train(int id, long long start_point, long long volume) {
  ifstream is(trainFile);
  is.seekg(start_point);

  string line;
  int sign = 1;
  int label;
  while (volume > 0 && getline(is, line)) {
    volume -= line.size() + 1;
    label = line[line.length() - 1] - '0';
    label_1_arry[id] += label;
    int count = 0, j = 0, sign = 1;
    if (label) {
      for (j = 0; count < FEATRURENS_NUM;) {
        if (line[j] == '-') {
          sign = -1;
          j++;
        }
        mean_fearture_1[id][count++] += (sign * (line[j + 2] - '0') * 100);
        j += 6;
        sign = 1;
      }
    } else {
      for (j = 0; count < FEATRURENS_NUM;) {
        if (line[j] == '-') {
          sign = -1;
          j++;
        }
        mean_fearture_0[id][count++] += (sign * (line[j + 2] - '0') * 100);
        j += 6;
        sign = 1;
      }
    }
    ++train_cnt[id];
  }
#ifdef TEST
  cout << "train:" << id << " time: "
       << duration_cast<duration<double>>(steady_clock::now() - start).count()
       << endl;
#endif
}

//找到每个线程读入数据的起始点
vector<pair<long long, long long>> get_train_point() {
  ifstream is(trainFile);
  long long length = FILE_SIZE;

  vector<pair<long long, long long>> res(THREADS_TRAIN);
  res[0].first = 0;
  string tmp;
  for (int i = 1; i < THREADS_TRAIN; ++i) {
    is.seekg(length / THREADS_TRAIN * i);
    getline(is, tmp);
    res[i].first = is.tellg();
    res[i - 1].second = res[i].first - res[i - 1].first;
  }

  is.seekg(length, ios::beg);
  getline(is, tmp);
  res.back().second = is.tellg();
  res.back().second -= (res.back().first + 1);

  is.close();

  // for (auto x : res)
  // 	printf("%lld, %lld %lld\n", x.first, x.second, x.first + x.second);

  return res;
}

vector<pair<long long, long long>> get_test_point() {
  ifstream is(testFile);
  is.seekg(0, is.end);
  long long length = is.tellg();

  vector<pair<long long, long long>> res(THREADS_TEST);

  res[0].first = 0;
  string tmp;

  for (int i = 1; i < THREADS_TEST; ++i) {
    is.seekg(length / THREADS_TEST * i);
    getline(is, tmp);
    res[i].first = is.tellg();
    res[i - 1].second = res[i].first - res[i - 1].first;
  }
  res.back().second = length - res.back().first;

  is.close();

  // for (auto x : res)
  // 	printf("%lld, %lld %lld\n", x.first, x.second, x.first + x.second);

  return res;
}

void StorePredictByString(char *res[THREADS_TEST]) {
  ofstream fout(predictFile);
  for (int i = 0; i < THREADS_TEST; i++) fout << res[i];
  fout.close();
}

int main(void) {
  ios::sync_with_stdio(false);
  vector<pair<long long, long long>> train_point = get_train_point();
  vector<pair<long long, long long>> test_point = get_test_point();
#ifdef TEST
  cout << "init read point time: "
       << duration_cast<duration<double>>(steady_clock::now() - start).count()
       << endl;
#endif
  vector<std::thread> test_thread;
  for (int i = 0; i < THREADS_TEST; ++i) {
    test_thread.push_back(
        thread(thread_read_test, i, test_point[i].first, test_point[i].second));
  }
  vector<std::thread> train_thread;
  for (int i = 0; i < THREADS_TRAIN; ++i) {
    train_thread.push_back(thread(thread_read_train, i, train_point[i].first,
                                  train_point[i].second));
  }

#ifdef TEST
  cout << "init thread time: "
       << duration_cast<duration<double>>(steady_clock::now() - start).count()
       << endl;
#endif
  for (auto &th : train_thread) {
    th.join();
  }
#ifdef TEST
  cout << "read train file time: "
       << duration_cast<duration<double>>(steady_clock::now() - start).count()
       << endl;
#endif

  int label_1 = accumulate(label_1_arry, label_1_arry + THREADS_TRAIN, 0);
  int line_num = accumulate(train_cnt, train_cnt + THREADS_TRAIN, 0);
  int label_0 = line_num - label_1;  // 2949 5051
  for (int i = 0; i < FEATRURENS_NUM; i++) {
    mean_fearture_1[0][i] +=
        mean_fearture_1[1][i] + mean_fearture_1[2][i] + mean_fearture_1[3][i];
    mean_fearture_0[0][i] +=
        mean_fearture_0[1][i] + mean_fearture_0[2][i] + mean_fearture_0[3][i];
    mean_fearture_1[0][i] /= label_1;
    mean_fearture_0[0][i] /= label_0;
  }
#ifdef TEST
  cout << "deal time: "
       << duration_cast<duration<double>>(steady_clock::now() - start).count()
       << endl;
#endif
  for (auto &th : test_thread) {
    th.join();
  }
#ifdef TEST
  cout << "read test time: "
       << duration_cast<duration<double>>(steady_clock::now() - start).count()
       << endl;
#endif
  char *res[THREADS_TEST];
  for (int i = 0; i < THREADS_TEST; i++) {
    res[i] = (char *)malloc(sizeof(char) * 2 * test_cnt[i] + 1);
    for (int j = 0; j < test_cnt[i]; j++) {
      res[i][j * 2] =
          get_distance(mean_fearture_1[0], mean_fearture_0[0], test_x[i][j]) +
          '0';
      res[i][j * 2 + 1] = '\n';
    }
    res[i][2 * test_cnt[i]] = '\0';
  }
  StorePredictByString(res);
#ifdef TEST
  cout << "all time: "
       << duration_cast<duration<double>>(steady_clock::now() - start).count()
       << endl;

  vector<int> pred;
  for (int i = 0; i < THREADS_TEST; i++) {
    for (int j = 0; j < test_cnt[i]; j++)
      pred.push_back(
          get_distance(mean_fearture_1[0], mean_fearture_0[0], test_x[i][j]));
  }
  vector<int> test_y;
  LoadTestLabels(test_y);
  cout << "acc: " << get_accuary(pred, test_y) << endl;
#endif
  return 0;
}
