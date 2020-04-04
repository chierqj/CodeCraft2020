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
#include <sstream>
#include <stack>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;

#ifdef TEST
#define TRAIN "../data/3512444/test_data.txt"
#define RESULT "../data/3512444/answer.txt"
#else
#define TRAIN "/data/test_data.txt"
#define RESULT "/projects/student/result.txt"
#endif

int del[4], mi = 0;
class FastFind {
 private:
  static const int MAXN = 560000 + 5;
  stack<int> sk;
  vector<bool> used;
  vector<int> beNum;
  vector<int> head;
  vector<int> preHead;
  vector<int> belong;
  vector<int> low;
  vector<int> dfn;
  vector<bool> vis;
  vector<int> reMap;
  vector<int> path;
  vector<int> degree;
  vector<vector<int>> answer[8];

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

 private:
  int edgeNum;
  int pointNum;
  int cnt;
  int SCCNum;
  int start;
  struct node {
    int to, nex, val;
  };
  int total = 0;
  unordered_map<int, int> mp;
  vector<bool> con[4];
  vector<node> edge;
  vector<node> pre;
};

void FastFind::Add(int x, int y, int val) {
  if (mp.find(x) == mp.end()) {
    mp[x] = ++pointNum;
  }
  if (mp.find(y) == mp.end()) {
    mp[y] = ++pointNum;
  }
  x = mp[x];
  y = mp[y];
  //正向边
  edge[++edgeNum].nex = head[x];
  edge[edgeNum].to = y;
  edge[edgeNum].val = val;
  head[x] = edgeNum;
  //反向边
  pre[edgeNum].nex = preHead[y];
  pre[edgeNum].to = x;
  pre[edgeNum].val = val;
  preHead[y] = edgeNum;
}
void FastFind::Init() {
  mp.clear();
  pointNum = 0;
  SCCNum = 0;
  cnt = 0;
  edgeNum = 0;
  head = vector<int>(MAXN, 0);
  preHead = vector<int>(MAXN, 0);
  reMap = vector<int>(MAXN);
  path = vector<int>(7);
  edge = vector<node>(MAXN);
  pre = vector<node>(MAXN);
  this->LoadData();
  vis = vector<bool>(pointNum + 1, false);
  con[1] = vector<bool>(pointNum + 1, false);
  con[2] = vector<bool>(pointNum + 1, false);
  con[3] = vector<bool>(pointNum + 1, false);
  low = vector<int>(pointNum + 1, 0);
  dfn = vector<int>(pointNum + 1, 0);
  beNum = vector<int>(pointNum + 1, 0);
  used = vector<bool>(pointNum + 1, false);
  degree = vector<int>(pointNum + 1, 0);
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
  if (degree[x] > 5) return true;
  if (con[1][x]) return true;
  for (int i = head[x]; i; i = edge[i].nex) {
    if (con[3][edge[i].to]) {
      return true;
    }
  }
  del[2]++;
  return false;
}
void FastFind::Tarjan(int x, int f) {
  low[x] = dfn[x] = ++cnt;
  vis[x] = true;
  sk.push(x);
  for (int i = head[x]; i; i = edge[i].nex) {
    int v = edge[i].to;
    if (v == f) continue;
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
  int fa;
  int tot = 0;
  for (int i = 1; i <= pointNum; i++) {
    while (belong[edge[head[i]].to] != belong[i]) {
      head[i] = edge[head[i]].nex;
      tot++;
      if (head[i] == 0) break;
    }
    fa = head[i];
    for (int j = edge[fa].nex; j; j = edge[j].nex) {
      if (belong[edge[j].to] == belong[i]) {
        edge[fa].nex = j;
        fa = j;
        degree[i]++;
      } else
        tot++;
    }
    edge[fa].nex = 0;
  }
  for (int i = 1; i <= pointNum; i++) {
    while (belong[pre[preHead[i]].to] != belong[i]) {
      preHead[i] = pre[preHead[i]].nex;
      tot++;
      if (preHead[i] == 0) break;
    }
    fa = preHead[i];
    for (int j = pre[fa].nex; j; j = pre[j].nex) {
      if (belong[pre[j].to] == belong[i]) {
        pre[fa].nex = j;
        fa = j;
      } else
        tot++;
    }
    pre[fa].nex = 0;
  }
  cout << tot << "\n";
}
void FastFind::FindCircle(int x, int dep) {
  path[dep] = x;
  int v;
  if (dep == 6) {
    if (con[1][x]) {
      answer[dep + 1].emplace_back(path);
      total++;
    }
    if (total % 10000 == 0) {
      cout << total << "\n";
    }
    return;
  }
  vis[x] = true;
  for (int i = head[x]; i; i = edge[i].nex) {
    v = edge[i].to;
    if (used[v]) {
      continue;
    }
    if (v == start && dep > 1) {
      answer[dep + 1].emplace_back(path);
      total++;
      if (total % 10000 == 0) {
        cout << total << "\n";
      }
      continue;
    }
    if (dep == 6) continue;
    if (vis[v]) continue;
    if (dep >= 3) {
      if (con[6 - dep][v]) this->FindCircle(v, dep + 1);

    } else {
      // if (dep == 2) {
      //     if (this->Judge(v)) this->FindCircle(v, dep + 1);
      // } else
      this->FindCircle(v, dep + 1);
    }
  }
  vis[x] = false;
}
void FastFind::Connect(int x, int dep) {
  if (dep == 3) return;
  vis[x] = true;
  for (int i = preHead[x]; i; i = pre[i].nex) {
    int v = pre[i].to;
    if (vis[v] || belong[v] != belong[x]) continue;
    for (int j = dep + 1; j <= 3; j++) con[j][v] = true;
    this->Connect(v, dep + 1);
  }
  vis[x] = false;
}
void FastFind::Search() {
  int x;
  cout << "-------------------\n";
  for (int i = 1; i <= pointNum; i++) {
    x = belong[i];
    if (beNum[x] > 2) {
      for (int j = 1; j <= 3; j++) {
        con[j] = vector<bool>(pointNum + 1, false);
      }
      Connect(i, 0);
      start = i;
      FindCircle(i, 0);
    }
    used[i] = true;
  }
}
void FastFind::WriteAnswer() {
  std::ofstream fout(RESULT);
  fout << total << "\n";
  for (int i = 2; i < 7; ++i) {
    for (auto &it : answer[i]) {
      for (int j = 0; j <= i; j++) {
        if (j != 0) fout << ",";
        fout << it[j];
      }
      fout << "\n";
    }
  }
  fout.close();
}
void FastFind::Run() {
  cout << "------开始初始化-------\n";
  cout << "---------------------\n";
  this->Init();
  cout << "---------------------\n";
  cout << "------初始化完成-------\n";
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
  cout << total << "\n";
  // WriteAnswer();
}
int main() {
  cout << "start\n";
  del[0] = 0;
  del[1] = 0;
  del[2] = 0;
  FastFind func;
  func.Run();
  // cout << mi << " " << del[2] << "\n";
  return 0;
}