#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/ipc.h>
#include <sys/mman.h>
#include <sys/shm.h>
#include <sys/stat.h>
#include <sys/wait.h>
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
#include <map>
#include <numeric>
#include <queue>
#include <sstream>
#include <stack>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;

#ifdef TEST
#define TRAIN "../data/std/test_data.txt"
#define RESULT "../data/std/result.txt"
#else
#define TRAIN "/data/test_data.txt"
#define RESULT "/projects/student/result.txt"
#endif

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
class FastFind {
 public:
  static const int MAXN = 560000 + 5;
  struct stringID {
    char id[10];
    int len = 0;
    int start;
  };

 public:
  void Run();
  void LoadData();
  void Init();
  void Tarjan(int x, int f);
  void Search();
  void WriteAnswer();
  void Add(int x, int y, int val);
  bool Judge(int &x);
  void GetSCC();
  void FindCircle(int x, int dep);
  void Connect(int x, int dep);
  void Delete();
  void Sort();
  void toString(int x, stringID &s);
  void ConnectBFS(int x);

 private:
  int edgeNum, pointNum;
  int cnt, SCCNum;
  int start;
  int total;
  struct node {
    int to, val;
  };
  stack<int> sk;
  unordered_map<int, int> mp;
  vector<bool> con[5];
  vector<node> G[MAXN];
  vector<node> Pre[MAXN];
  stringID outID[MAXN];
  int ID[MAXN], dfsID[MAXN];
  vector<bool> used;
  vector<int> beNum, belong, low, dfn;
  vector<bool> vis;
  vector<int> path;
  vector<vector<int>> answer[8];
  map<int, int> fk;
  int mxNum = 0;
  int del = 0;
  vector<int> newG[MAXN];
  vector<bool> pointVis;
  vector<bool> conVis;
};

void FastFind::Add(int x, int y, int val) {
  if (mp.find(x) == mp.end()) {
    mp[x] = ++pointNum;
    ID[pointNum] = x;
  }
  if (mp.find(y) == mp.end()) {
    mp[y] = ++pointNum;
    ID[pointNum] = y;
  }
  x = mp[x];
  y = mp[y];
  edgeNum++;
  //正向边
  G[x].emplace_back(node{y, val});
  //反向边
  Pre[y].emplace_back(node{x, val});
}
void FastFind::Init() {
  mp.clear();
  pointNum = 0;
  SCCNum = 0;
  total = 0;
  cnt = 0;
  edgeNum = 0;
  path = vector<int>(7);
  this->LoadData();
  vis = vector<bool>(pointNum + 1, false);
  con[1] = vector<bool>(pointNum + 1, false);
  con[2] = vector<bool>(pointNum + 1, false);
  con[3] = vector<bool>(pointNum + 1, false);
  low = vector<int>(pointNum + 1, 0);
  dfn = vector<int>(pointNum + 1, 0);
  beNum = vector<int>(pointNum + 1, 0);
  used = vector<bool>(pointNum + 1, false);
  belong = vector<int>(pointNum + 1);
}
void FastFind::LoadData() {
  struct stat sb;
  int fd = open(TRAIN, O_RDONLY);
  fstat(fd, &sb);
  long long bufferSize = sb.st_size;
  char *buffer = (char *)mmap(NULL, bufferSize, PROT_READ, MAP_PRIVATE, fd, 0);

  close(fd);
  int num = 0, points[3], index = 0;
  long long id = 0;
  while (id < bufferSize) {
    if (buffer[id] == ',') {
      points[index++] = num;
      num = 0;
    } else if (buffer[id] == '\n') {
      points[index++] = num;
      this->Add(points[0], points[1], points[2]);
      index = 0;
      num = 0;
    } else {
      num = num * 10 + *(buffer + id) - '0';
    }
    id++;
  }
}

bool FastFind::Judge(int &x) {
  if (con[1][x]) return true;
  for (auto &it : G[x]) {
    if (con[3][it.to]) {
      return true;
    }
  }
  return false;
}
void FastFind::Tarjan(int x, int f) {
  low[x] = dfn[x] = ++cnt;
  vis[x] = true;
  sk.push(x);
  for (auto &it : G[x]) {
    int v = it.to;
    // if (v == f) continue;
    if (!dfn[v]) {
      Tarjan(v, x);
      low[x] = min(low[x], low[v]);
    } else if (vis[v]) {
      low[x] = min(low[x], dfn[v]);
    }
  }
  if (low[x] == dfn[x]) {
    SCCNum++;
    while (1) {
      int p = sk.top();
      sk.pop();
      belong[p] = SCCNum;
      vis[p] = false;
      beNum[SCCNum]++;
      if (p == x) break;
    }
  }
}
void FastFind::GetSCC() {
  int n = pointNum;
  for (int i = 1; i <= n; i++) {
    if (!dfn[i]) this->Tarjan(i, 0);
  }
}

void FastFind::Delete() {
  for (int i = 1; i <= pointNum; i++) {
    int v = belong[i];
    for (auto it = G[i].begin(); it != G[i].end();) {
      if (belong[(*it).to] != v) {
        del++;
        it = G[i].erase(it);
      } else {
        ++it;
      }
    }
    for (auto it = Pre[i].begin(); it != Pre[i].end();) {
      if (belong[(*it).to] != v) {
        it = Pre[i].erase(it);
      } else {
        ++it;
      }
    }
  }
}
void FastFind::FindCircle(int x, int dep) {
  path[dep] = x;
  int v;
  if (dep == 6) {
    //走一步可以到start
    if (con[1][x]) {
      answer[dep + 1].emplace_back(path);
    }
    return;
  }
  vis[x] = true;
  if (dep > 3) {
    if (conVis[x]) {
      sort(newG[x].begin(), newG[x].end(), [&](const int &a, const int &b) {
        return this->ID[a] < this->ID[b];
      });
      conVis[x] = false;
    }
    for (auto &it : newG[x]) {
      if (used[it]) {
        continue;
      }
      if (it == start) {
        answer[dep + 1].emplace_back(path);
        continue;
      }
      if (vis[it]) continue;
      if (con[6 - dep][it]) this->FindCircle(it, dep + 1);
    }
  } else {
    for (auto &it : G[x]) {
      v = it.to;
      if (used[v]) {
        continue;
      }
      if (v == start && dep > 1) {
        answer[dep + 1].emplace_back(path);
        continue;
      }
      if (vis[v]) continue;
      if (dep > 2) {
        if (con[6 - dep][v]) this->FindCircle(v, dep + 1);
      } else {
        if (dep == 2) {
          if (this->Judge(v)) this->FindCircle(v, dep + 1);
        } else
          this->FindCircle(v, dep + 1);
      }
    }
  }
  vis[x] = false;
}

void FastFind::Connect(int x, int dep) {
  if (dep == 3) return;
  vis[x] = true;
  for (auto &it : Pre[x]) {
    int v = it.to;
    if (!conVis[x]) {
      if (pointVis[v])
        newG[v].emplace_back(x);
      else {
        pointVis[v] = true;
        newG[v].clear();
        newG[v].emplace_back(x);
      }
    }
    if (vis[v]) continue;
    for (int j = dep + 1; j <= 3; j++) con[j][v] = true;
    this->Connect(v, dep + 1);
  }
  conVis[x] = true;
  vis[x] = false;
}
void FastFind::Sort() {
  for (int i = 1; i <= pointNum; i++) dfsID[i] = i;
  sort(dfsID + 1, dfsID + 1 + pointNum,
       [&](const int &a, const int &b) { return this->ID[a] < this->ID[b]; });
  for (int i = 1; i <= pointNum; i++) {
    sort(G[i].begin(), G[i].end(), [&](const node &a, const node &b) {
      return this->ID[a.to] < this->ID[b.to];
    });
  }
}
void FastFind::Search() {
  int x;
  this->Sort();
  for (int i = 1; i <= pointNum; i++) {
    x = belong[dfsID[i]];
    if (beNum[x] > 2) {
      conVis = vector<bool>(pointNum + 1, false);
      pointVis = vector<bool>(pointNum + 1, false);
      for (int j = 1; j <= 3; j++) {
        con[j] = vector<bool>(pointNum + 1, false);
      }
      Connect(dfsID[i], 0);
      conVis = vector<bool>(pointNum + 1, true);
      start = dfsID[i];
      FindCircle(dfsID[i], 0);
    }
    used[dfsID[i]] = true;
  }
}

void FastFind::toString(int x, stringID &s) {
  if (x == 0) {
    s.id[9] = '0';
    s.len = 1;
    s.start = 9;
    return;
  }
  s.start = 10;
  while (x) {
    s.id[--s.start] = x % 10 + '0';
    x /= 10;
  }
  s.len = 10 - s.start;
}

void FastFind::WriteAnswer() {
  for (int i = 1; i <= pointNum; i++) {
    this->toString(ID[i], outID[i]);
  }
  uint32_t bufsize = 80 * total + 15;
  uint32_t idx = 0;
  int v;
  char *buffer = new char[bufsize];
  stringID tot;
  this->toString(total, tot);
  memcpy(buffer + idx, tot.id + tot.start, tot.len);
  idx += tot.len;
  buffer[idx++] = '\n';
  for (int len = 3; len < 8; ++len) {
    for (int i = 0; i < answer[len].size(); ++i) {
      const auto &row = answer[len][i];
      for (int k = 0; k < len; ++k) {
        v = row[k];
        if (k != 0) buffer[idx++] = ',';
        memcpy(buffer + idx, outID[v].id + outID[v].start, outID[v].len);
        idx += outID[v].len;
      }
      buffer[idx++] = '\n';
    }
  }
  int fd = open(RESULT, O_RDWR | O_CREAT, 0666);
  char *result =
      (char *)mmap(NULL, idx, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
  ftruncate(fd, idx);
  close(fd);
  memcpy(result, buffer, idx);
}
void FastFind::Run() {
  // cout << "------开始初始化-------\n";
  // cout << "---------------------\n";
  this->Init();
  // cout << "---------------------\n";
  // cout << "------初始化完成-------\n";
  // cout << "--节点数量: " << pointNum << " --\n";
  // cout << "--边数量: " << edgeNum << " --\n";
  // cout << "---------------------\n";
  // cout << "------切分强连通分量-------\n";
  // cout << "---------------------\n";
  this->GetSCC();
  // cout << "---------------------\n";
  // cout << "------切分完成-------\n";
  // cout << "--"
  //      << "强连通个数: " << SCCNum << " --\n";
  this->Delete();
  // cout << "---------------------\n";
  // cout << "------搜索环-------\n";
  // cout << "---------------------\n";
  this->Search();
  for (int i = 3; i <= 7; i++) {
    total += answer[i].size();
    cout << answer[i].size() << "\n";
  }
  cout << "----del: " << del << " ---\n";
  cout << "----环个数: " << total << " -----\n";
  ScopeTime t;
  WriteAnswer();
  t.LogTime();
}

int main() {
  FastFind *func = new FastFind();
  func->Run();
  return 0;
}