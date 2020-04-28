#include <chrono>
#include <iostream>
#include <random>
#include <set>
#include <vector>
using namespace std;
int main() {
  srand((unsigned long long)new char);
  set<pair<int, int>> vis;
  vector<pair<int, int>> E;
  auto add_edge = [&](int x, int y) -> bool {
    if (x == y || vis.count(make_pair(x, y))) return false;
    vis.insert(make_pair(x, y));
    E.push_back(make_pair(x, y));
    return true;
  };
  int e = 280000;
  for (int i = 50001; i <= 200000; i++) {
    add_edge(rand() % (i - 50000) + 50000, i);
    e -= 1;
  }
  for (int i = 6000; i <= 6000 + 11; i++) {
    for (int j = 6000; j <= 6000 + 11; j++) {
      if (i == j) continue;
      add_edge(i, j);
      e -= 1;
    }
  }
  for (int i = 10000; i <= 10000 + 12; i++) {
    for (int j = 10000; j <= 10000 + 12; j++) {
      if (i == j) continue;
      add_edge(i, j);
      e -= 1;
    }
  }
  for (int i = 25123; i <= 25123 + 11; i++) {
    for (int j = 25123; j <= 25123 + 11; j++) {
      if (i == j) continue;
      add_edge(i, j);
      e -= 1;
    }
  }
  while (e) {
    int t = rand() % 450 + 5000;
    int x = rand() % (50000 - t);
    int y = x + t;
    if (rand() & 1) swap(x, y);
    if (add_edge(x, y)) e -= 1;
  }
  random_shuffle(E.begin(), E.end());
  for (auto e : E) {
    cout << e.first << "," << e.second << ",123123\n";
  }
  return 0;
}