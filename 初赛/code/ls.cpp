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
#define TRAIN "../data/1004812/test_data.txt"
#define RESULT "../data/1004812/result.txt"
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
  static const int THREADNUM = 4;
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
  bool Judge(int &x, int &pid);
  void GetSCC();
  void FindCircle(int x, int dep, int &pid);
  void Connect(int x, int dep, int &pid);
  void Delete();
  void Sort();
  void toString(int x, stringID &s);
  void ConnectBFS(int x);
  void ThreadSearch(int startId, int endId, int pid);
  void ThreadToChar(int pid);
  void ThreadRun(int startId, int endId, int pid);

 private:
  int edgeNum, pointNum;
  int cnt, SCCNum;
  int total[THREADNUM];
  struct node {
    int to, val;
  };
  stack<int> sk;
  unordered_map<int, int> mp;
  vector<node> G[MAXN];
  vector<node> Pre[MAXN];
  stringID outID[MAXN];
  int ID[MAXN], dfsID[MAXN];
  int del = 0;

 public:
  vector<bool> used[THREADNUM];
  vector<int> beNum, belong, low, dfn;
  vector<bool> vis[THREADNUM];
  vector<int> path[THREADNUM];
  vector<vector<int>> answer[THREADNUM][8];
  int idx[THREADNUM][8];
  char *buf[THREADNUM][8];
  vector<bool> con[THREADNUM][4];
  int start[THREADNUM];
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
  cnt = 0;
  edgeNum = 0;
  this->LoadData();
  vis[0] = vector<bool>(pointNum + 1, false);
  for (int i = 0; i < THREADNUM; i++) {
    path[i] = vector<int>(7);
    con[i][1] = vector<bool>(pointNum + 1, false);
    con[i][2] = vector<bool>(pointNum + 1, false);
    con[i][3] = vector<bool>(pointNum + 1, false);
    used[i] = vector<bool>(pointNum + 1, false);
  }
  low = vector<int>(pointNum + 1, 0);
  dfn = vector<int>(pointNum + 1, 0);
  beNum = vector<int>(pointNum + 1, 0);
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

bool FastFind::Judge(int &x, int &pid) {
  if (con[pid][1][x]) return true;
  for (auto &it : G[x]) {
    if (con[pid][3][it.to]) {
      return true;
    }
  }
  return false;
}
void FastFind::Tarjan(int x, int f) {
  low[x] = dfn[x] = ++cnt;
  vis[0][x] = true;
  sk.push(x);
  for (auto &it : G[x]) {
    int v = it.to;
    // if (v == f) continue;
    if (!dfn[v]) {
      Tarjan(v, x);
      low[x] = min(low[x], low[v]);
    } else if (vis[0][v]) {
      low[x] = min(low[x], dfn[v]);
    }
  }
  if (low[x] == dfn[x]) {
    SCCNum++;
    while (1) {
      int p = sk.top();
      sk.pop();
      belong[p] = SCCNum;
      vis[0][p] = false;
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
void FastFind::FindCircle(int x, int dep, int &pid) {
  path[pid][dep] = x;
  int v;
  if (dep == 6) {
    //走一步可以到start
    if (con[pid][1][x]) {
      answer[pid][dep + 1].emplace_back(path[pid]);
    }
    return;
  }

  vis[pid][x] = true;
  for (auto &it : G[x]) {
    v = it.to;
    if (used[pid][v]) {
      continue;
    }
    if (v == start[pid] && dep > 1) {
      answer[pid][dep + 1].emplace_back(path[pid]);
      continue;
    }
    if (vis[pid][v]) continue;
    if (dep > 2) {
      if (con[pid][6 - dep][v]) this->FindCircle(v, dep + 1, pid);
    } else {
      if (dep == 2) {
        if (this->Judge(v, pid)) this->FindCircle(v, dep + 1, pid);
      } else
        this->FindCircle(v, dep + 1, pid);
    }
  }
  vis[pid][x] = false;
}

void FastFind::Connect(int x, int dep, int &pid) {
  if (dep == 3) return;
  vis[pid][x] = true;
  for (auto &it : Pre[x]) {
    int v = it.to;
    if (vis[pid][v]) continue;
    for (int j = dep + 1; j <= 3; j++) con[pid][j][v] = true;
    this->Connect(v, dep + 1, pid);
  }
  vis[pid][x] = false;
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
void FastFind::ThreadSearch(int startId, int endId, int pid) {
  int x;
  vis[pid] = vector<bool>(pointNum + 1, false);
  for (int i = startId; i < endId; i++) {
    x = belong[dfsID[i]];
    if (beNum[x] > 2) {
      for (int j = 1; j <= 3; j++) {
        con[pid][j] = vector<bool>(pointNum + 1, false);
      }
      start[pid] = dfsID[i];
      Connect(dfsID[i], 0, pid);
      FindCircle(dfsID[i], 0, pid);
    }
    used[pid][dfsID[i]] = true;
  }
}
void FastFind::ThreadToChar(int pid) {
  uint32_t bufsize = 80 * total[pid] + 15;
  uint32_t id;
  int v;
  for (int i = 3; i < 8; i++) {
    buf[pid][i] = new char[bufsize];
  }
  for (int len = 3; len < 8; ++len) {
    id = 0;
    for (int i = 0; i < answer[pid][len].size(); ++i) {
      const auto &row = answer[pid][len][i];
      for (int k = 0; k < len; ++k) {
        v = row[k];
        if (k != 0) buf[pid][len][id++] = ',';
        memcpy(buf[pid][len] + id, outID[v].id + outID[v].start, outID[v].len);
        id += outID[v].len;
      }
      buf[pid][len][id++] = '\n';
    }
    idx[pid][len] = id;
  }
}
void FastFind::ThreadRun(int startId, int endId, int pid) {
  for (int i = 1; i < startId; i++) {
    used[pid][dfsID[i]] = true;
  }
  this->ThreadSearch(startId, endId, pid);
  total[pid] = 0;
  for (int i = 3; i <= 7; i++) {
    total[pid] += answer[pid][i].size();
  }
  this->ThreadToChar(pid);
}

void FastFind::Search() {
  int x;
  this->Sort();
  for (int i = 1; i <= pointNum; i++) {
    this->toString(ID[i], outID[i]);
  }
  int p[4] = {2, 3, 5, 40};
  vector<thread> Threads(THREADNUM);
  int pre = 1, top, block = (pointNum + 10) / 50;
  for (int i = 0; i < THREADNUM; i++) {
    top = pre + block * p[i];
    if (top > pointNum + 1) top = pointNum + 1;
    Threads[i] = thread(&FastFind::ThreadRun, this, pre, top, i);
    pre = top;
  }
  for (auto &it : Threads) it.join();
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
  int fd = open(RESULT, O_RDWR | O_CREAT, 0666);
  for (int i = 1; i < THREADNUM; i++) {
    total[0] += total[i];
  }
  uint32_t bufsize = 80 * total[0] + 15;
  uint32_t id = 0;
  char *buffer = new char[bufsize];
  stringID ansNum;
  toString(total[0], ansNum);
  memcpy(buffer, ansNum.id + ansNum.start, ansNum.len);
  id += ansNum.len;
  buffer[id] = '\n';
  id++;
  for (int i = 3; i < 8; i++) {
    for (int j = 0; j < THREADNUM; j++) {
      memcpy(buffer + id, buf[j][i], idx[j][i]);
      id += idx[j][i];
    }
  }
  char *result =
      (char *)mmap(NULL, id, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
  ftruncate(fd, id);
  close(fd);
  memcpy(result, buffer, id);
}
void FastFind::Run() {
  cout << "------开始初始化-------\n";
  cout << "---------------------\n";
  this->Init();
  cout << "---------------------\n";
  cout << "------初始化完成-------\n";
  cout << "--节点数量: " << pointNum << " --\n";
  cout << "--边数量: " << edgeNum << " --\n";
  cout << "---------------------\n";
  cout << "------切分强连通分量-------\n";
  cout << "---------------------\n";
  this->GetSCC();
  cout << "---------------------\n";
  cout << "------切分完成-------\n";
  cout << "--"
       << "强连通个数: " << SCCNum << " --\n";
  this->Delete();
  cout << "---------------------\n";
  cout << "------搜索环-------\n";
  cout << "---------------------\n";
  this->Search();
  int nu = 0;
  for (int i = 3; i < 8; i++) {
    nu = 0;
    for (int j = 0; j < THREADNUM; j++) {
      nu += answer[j][i].size();
    }
    cout << nu << "\n";
  }
  WriteAnswer();
  cout << "-------------\n";
  cout << total[0] << "\n";
}

int main() {
  FastFind *func = new FastFind();
  func->Run();
  return 0;
}