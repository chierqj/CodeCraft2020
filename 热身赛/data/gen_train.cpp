#include <iostream>
#include <random>
#include <vector>
using namespace std;

int main() {
  mt19937 rng(random_device{}());
  const int W = 80000;
  const int V = 1000;
  // const int W = 10;
  // const int V = 5;
  vector<double> cof(V);
  auto rnd = [&]() -> double { return (double)(rng() % 1001) / 1000; };
  for (int i = 0; i < V; i++) {
    cof[i] = rnd() * 2 - 1;
  }
  vector<pair<double, vector<double>>> vec;
  for (int i = 0; i < W; i++) {
    vector<double> cur(V);
    bool flag = false;
    for (int t = 0; t < V; t++) {
      if (flag) {
        cur[t] = rnd() * 2 - 1;
      } else {
        cur[t] = rnd();
      }
    }
    double ret = 0.0;
    for (int t = 0; t < V; t++) {
      ret += cur[t] * cof[t];
    }
    vec.push_back(make_pair(ret, cur));
  }
  sort(vec.begin(), vec.end());
  double magic = vec[(int)vec.size() / 3].first;
  random_shuffle(vec.begin(), vec.end());
  FILE* f1 = fopen("train_data.txt", "w");
  for (int i = 0; i < (int)vec.size(); i++) {
    auto& cur = vec[i].second;
    for (auto v : cur) {
      fprintf(f1, "%.3f,", v);
    }
    if (vec[i].first <= magic) {
      fprintf(f1, "0\n");
    } else {
      fprintf(f1, "1\n");
    }
  }
  fclose(f1);
  FILE* f2 = fopen("ddd.txt", "w");
  FILE* f3 = fopen("ddd_answer.txt", "w");
  // const int W2 = 10;
  // const int V2 = 5;
  const int W2 = 20000;
  const int V2 = 1000;
  cout << "!" << endl;
  for (int i = 0; i < W2; i++) {
    double sum = 0;
    for (int t = 0; t < V2; t++) {
      double v = rnd();
      sum += cof[t] * v;
      fprintf(f2, "%.3f", v);
      if (t + 1 == V2) {
        fprintf(f2, "\n");
      } else {
        fprintf(f2, ",");
      }
    }
    if (sum <= magic) {
      fprintf(f3, "0\n");
    } else {
      fprintf(f3, "1\n");
    }
  }
  fclose(f2);
  fclose(f3);
  return 0;
}
