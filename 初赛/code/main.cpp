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
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef LOCAL
#define TRAIN "../data/3512444/test_data.txt"
#define RESULT "../data/3512444/result.txt"
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

std::ostream &operator<<(std::ostream &os, const std::vector<int> &mat) {
  os << "{";
  for (int i = 0; i < mat.size(); ++i) {
    if (i != 0) os << ",";
    os << mat[i];
  }
  os << "}";
  return os;
}
std::ostream &operator<<(std::ostream &os,
                         const std::vector<std::vector<int>> &mat) {
  os << "[";
  for (int i = 0; i < mat.size(); ++i) {
    if (i != 0) os << ",";
    os << mat[i];
  }
  os << "]";
  return os;
}

/*
	本端账号ID和对端账号ID为一个32位的正整数
	转账金额为一个32位的正整数
	转账记录最多为28万条
	每个账号平均转账记录数< 10
	账号A给账号B最多转账一次
*/

class XJBG {
 public:
  static const int NTHREAD = 4;        // 线程个数
  static const int LIMIT_STEP = 3;     // 预估步长
  static const int MAXN = 560000 + 7;  // 总点数
  const int PARAM[NTHREAD] = {2, 3, 5, 7};

 private:
  std::unordered_map<int, int> m_IDToMap;  // ID->MapID
  std::vector<int> m_IDDom;                // MapID->ID
  int m_answers = 0;                       // 总环数目
  std::vector<int> m_Circles;              // 大于3的联通分量的点

  struct Node {
    int answers;                              // 环个数
    std::vector<int> tempPath;                // dfs临时路径
    std::vector<bool> stepArrive[3];          // 预估步长
    std::vector<std::vector<int>> Answer[5];  // 结果
    bool vis[MAXN];                           // 标记访问
    Node() { tempPath = std::vector<int>(7, -1); }
  };
  std::vector<Node> ThreadData;  // 线程数据

 private:
  struct Edge {
    int v, w;
  };
  std::vector<Edge> m_Sons[MAXN];               // 边
  std::vector<Edge> m_Fathers[MAXN];            // 边
  int m_edgeNum = 0;                            // 边数目
  int m_dfn[MAXN], m_low[MAXN], m_stack[MAXN];  // Tarjan
  int m_category[MAXN];                         // 所在联通分量id
  int m_inDegree[MAXN], m_outDegree[MAXN];      // 出入度
  bool m_inStack[MAXN];                         // Tarjan标记在不在栈内
  bool m_inCircle[MAXN];                        // 在不在找过的环内
  int m_tarjanCount = 0;                        // Tarjan搜索顺序
  int m_stackTop = 0;                           // Tarjan栈top标号
  int m_scc = 0, m_useScc = 0;                  // 总联通分量和大于3的

 public:
  void Init();
  void LoadData();
  void TarJan();
  void FindPath();
  void SaveAnswer();

 private:
  inline int GetID(int x);
  int GetMapID(int x);
  inline void addEdge(int u, int v, int w);
  void tarjan(int u);
  void pretreatTravelBFS(Node &Data, int &st);
  bool judge(Node &Data, const int &ctg, const int &u, const int &dep);
  void handelCircle(Node &Data, const int &dep);
  void findCircle(Node &Data, const int &ctg, int u, int dep);
};

void XJBG::Init() { ThreadData = std::vector<Node>(NTHREAD); }
inline int XJBG::GetID(int x) { return m_IDDom[x]; }
int XJBG::GetMapID(int x) {
  auto it = m_IDToMap.find(x);
  if (it != m_IDToMap.end()) {
    return it->second;
  }
  int sz = m_IDDom.size();
  m_IDToMap.insert({x, sz});
  m_IDDom.emplace_back(x);
  return sz;
}
inline void XJBG::addEdge(int u, int v, int w) {
  u = this->GetMapID(u);
  v = this->GetMapID(v);
  m_Sons[u].emplace_back(Edge{v, w});
  m_Fathers[v].emplace_back(Edge{u, w});
  ++m_edgeNum;
  ++m_inDegree[v];
  ++m_outDegree[u];
}

void XJBG::LoadData() {
  ScopeTime t;

  struct stat sb;
  int fd = open(TRAIN, O_RDONLY);
  fstat(fd, &sb);
  long long bufsize = sb.st_size;
  char *buffer = (char *)mmap(NULL, bufsize, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);

  int temp[2], idx = 0, x = 0;
  long long start = 0;
  while (start < bufsize) {
    if (*(buffer + start) == ',') {
      temp[idx++] = x;
      x = 0;
    } else if (*(buffer + start) == '\n') {
      addEdge(temp[0], temp[1], x);
      idx = 0;
      x = 0;
    } else {
      x = x * 10 + (*(buffer + start) - '0');
    }
    ++start;
  }

#ifdef LOCAL
  std::cerr << "@ LoadData: (V: " << m_IDDom.size() << ", E: " << m_edgeNum
            << ") #";
  t.LogTime();
#endif
#ifdef TEST
  std::set<std::vector<int>> st;
  for (int i = 0; i < MAXN; ++i) {
    for (auto &it : m_Sons[i]) {
      st.insert({i, it.v});
    }
  }
  std::cerr << "@ unique edge: " << st.size() << "\n";
#endif
}

void XJBG::tarjan(int u) {
  m_dfn[u] = m_low[u] = ++m_tarjanCount;
  m_stack[m_stackTop++] = u;
  m_inStack[u] = 1;

  for (auto &it : m_Sons[u]) {
    int v = it.v;
    if (!m_dfn[v]) {
      this->tarjan(v);
      m_low[u] = std::min(m_low[u], m_low[v]);
    } else if (m_inStack[v]) {
      m_low[u] = std::min(m_low[u], m_dfn[v]);
    }
  }

  if (m_dfn[u] == m_low[u]) {
    std::vector<int> tmp;
    int cur;
    do {
      cur = m_stack[--m_stackTop];
      m_inStack[cur] = 0;
      m_category[cur] = u;
      tmp.emplace_back(cur);
    } while (cur != u);
    if (tmp.size() >= 3) {
      ++m_useScc;
      m_Circles.insert(m_Circles.end(), tmp.begin(), tmp.end());
    }
    ++m_scc;
  }
}
void XJBG::TarJan() {
  ScopeTime t;

  for (int i = 0; i < m_IDDom.size(); ++i) {
    if (!m_dfn[i] && !m_Sons[i].empty()) {
      this->tarjan(i);
    }
  }
  std::sort(m_Circles.begin(), m_Circles.end(),
            [&](const int &x, const int &y) {
              return this->GetID(x) < this->GetID(y);
            });

#ifdef LOCAL
  std::cerr << "@ TarJan: (scc: " << m_scc << ", usescc: " << m_useScc << ") #";
  t.LogTime();
#endif
}

void XJBG::pretreatTravelBFS(Node &Data, int &st) {
  std::queue<std::pair<int, int>> Q;
  std::vector<bool> vis(m_IDDom.size(), false);
  Q.push(std::make_pair(st, 0));
  vis[st] = true;
  int ctg = m_category[st];
  while (!Q.empty()) {
    const auto &head = Q.front();
    Q.pop();
    if (head.second == 1) {
      Data.stepArrive[0][head.first] = true;
      Data.stepArrive[1][head.first] = true;
      Data.stepArrive[2][head.first] = true;
    } else if (head.second == 2) {
      Data.stepArrive[1][head.first] = true;
      Data.stepArrive[2][head.first] = true;
    } else if (head.second == 3) {
      Data.stepArrive[2][head.first] = true;
      continue;
    }
    for (auto &it : m_Fathers[head.first]) {
      int v = it.v;
      if (vis[v]) continue;
      if (m_category[v] != ctg) {
        Data.vis[v] = true;
        continue;
      }
      vis[v] = true;
      Q.push(std::make_pair(v, head.second + 1));
    }
  }
}

void XJBG::handelCircle(Node &Data, const int &dep) {
  Data.Answer[dep - 3].emplace_back(Data.tempPath);
  ++Data.answers;
}
bool XJBG::judge(Node &Data, const int &ctg, const int &v, const int &dep) {
  if (Data.vis[v] || dep > 6) return true;
  if (dep == 4 && !Data.stepArrive[2][v]) {
    return true;
  }
  if (dep == 5 && !Data.stepArrive[1][v]) {
    return true;
  }
  if (dep == 6) {
    if (Data.stepArrive[0][v]) this->handelCircle(Data, dep + 1);
    return true;
  }
  return false;
}
void XJBG::findCircle(Node &Data, const int &ctg, int u, int dep) {
  for (auto &it : m_Sons[u]) {
    int v = it.v;
    Data.tempPath[dep] = v;
    if (v == Data.tempPath[0] && dep >= 3) {
      this->handelCircle(Data, dep);
      continue;
    }
    if (this->judge(Data, ctg, v, dep)) {
      continue;
    }
    Data.vis[v] = true;
    this->findCircle(Data, ctg, v, dep + 1);
    Data.vis[v] = false;
  }
}

void XJBG::FindPath() {
  ScopeTime t;

  for (auto &v : m_Circles) {
    std::sort(m_Sons[v].begin(), m_Sons[v].end(),
              [&](const Edge &e1, const Edge &e2) {
                return this->GetID(e1.v) < this->GetID(e2.v);
              });
  }

  auto foo = [&](int pid, int start, int end) {
    auto &Data = ThreadData[pid];
    for (int i = 0; i < start; ++i) Data.vis[m_Circles[i]] = true;
    for (int i = start; i < end; ++i) {
      for (auto &it : Data.stepArrive) {
        it = std::vector<bool>(m_IDDom.size(), false);
      }
      int v = m_Circles[i];
      this->pretreatTravelBFS(Data, v);
      Data.tempPath[0] = v;
      Data.vis[v] = true;
      this->findCircle(Data, m_category[v], v, 1);
    }
  };

  int sum = 0;
  for (int i = 0; i < NTHREAD; ++i) sum += PARAM[i];
  int start = 0, block = m_Circles.size() / sum;
  std::vector<std::thread> Threads(NTHREAD);
  for (int i = 0; i < NTHREAD; ++i) {
    int end = (i == NTHREAD - 1 ? m_Circles.size() : start + block * PARAM[i]);
    Threads[i] = std::thread(foo, i, start, end);
    start += block * PARAM[i];
  }
  for (auto &it : Threads) it.join();
  for (auto &it : ThreadData) m_answers += it.answers;

#ifdef LOCAL
  std::cerr << "@ Answer: " << m_answers << " #";
  t.LogTime();
  for (int i = 0; i < NTHREAD; ++i) {
    const auto &data = ThreadData[i];
    std::cerr << "* " << i << ": " << data.answers << "(";
    for (int j = 0; j < 4; ++j) {
      std::cerr << data.Answer[j].size() << ", ";
    }
    std::cerr << data.Answer[4].size() << ")\n";
  }
#endif
}

void XJBG::SaveAnswer() {
  struct T {
    char str[10];
    int length = 0;
  };
  std::vector<T> MapID(m_IDDom.size());
  for (int i = 0; i < m_IDDom.size(); ++i) {
    int x = m_IDDom[i];
    if (x == 0) {
      MapID[i].length = 1;
      MapID[i].str[0] = '0';
      continue;
    }
    int idx = 0;
    while (x) {
      MapID[i].str[idx++] = x % 10 + '0';
      x /= 10;
    }
    MapID[i].length = idx;
  }

  struct T1 {
    char *buffer[5];
    int length[5];
  };
  std::vector<T1> OutData(NTHREAD);

  // auto foo = [&](int pid) {
  //   const auto &data = ThreadData[pid];
  //   auto &out = OutData[pid];
  //   for (int len = 4; len >= 0; --len) {
  //     out.buffer[len] = new char[data.Answer[len].size() * 80];
  //     int idx = 0;
  //     for (int i = data.Answer[len].size() - 1; i >= 0; --i) {
  //       const auto &row = data.Answer[len][i];
  //       buffer[idx++] = '\n';
  //       for (int k = len + 2; k >= 0; --k) {
  //         for (int kk = 0; kk < MapID[row[k]].length; ++kk) {
  //           buffer[idx++] = MapID[row[k]].str[kk];
  //         }
  //         if (k != 0) buffer[idx++] = ',';
  //       }
  //     }
  //   }
  // };

  ScopeTime t;

  uint32_t bufsize = 80 * m_answers + 15;
  uint32_t idx = 0;
  char *buffer = new char[bufsize];

  int x = m_answers;
  if (x == 0) {
    buffer[idx++] = '0';
  } else {
    char tmp[10];
    int ti = 0;
    while (x) {
      tmp[ti++] = x % 10 + '0';
      x /= 10;
    }
    for (int i = ti - 1; i >= 0; --i) buffer[idx++] = tmp[i];
  }
  buffer[idx++] = '\n';
  for (int len = 0; len < 5; ++len) {
    for (auto &data : ThreadData) {
      const auto &ans = data.Answer[len];
      for (auto &it : ans) {
        // 便利一行len+3个数字
        for (int i = 0; i < len + 3; ++i) {
          const int &num = it[i];
          if (i != 0) buffer[idx++] = ',';
          for (int j = MapID[num].length - 1; j >= 0; --j) {
            buffer[idx++] = MapID[num].str[j];
          }
        }
        buffer[idx++] = '\n';
      }
    }
  }

  int fd = open(RESULT, O_RDWR | O_CREAT, 0666);
  char *result =
      (char *)mmap(NULL, idx, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
  ftruncate(fd, idx);
  close(fd);
  memcpy(result, buffer, idx);

#ifdef LOCAL
  std::cerr << "@ WriteAnswer # ";
  t.LogTime();
#endif
}

int main() {
  std::cerr << std::fixed << std::setprecision(3);

  XJBG *xjbg = new XJBG();
  xjbg->Init();
  xjbg->LoadData();
  xjbg->TarJan();
  xjbg->FindPath();
  xjbg->SaveAnswer();
  return 0;
}