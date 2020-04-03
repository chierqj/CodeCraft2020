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
 private:
  static const int MAXN = 560000 + 7;      // 总点数
  std::unordered_map<int, int> m_IDToMap;  // ID->MapID
  std::vector<int> m_IDDom;                // MapID->ID

  int m_answers = 0;                          // 总环数目
  std::vector<int> m_Circles;                 // 大于3的联通分量的点
  std::vector<int> m_TempPath;                // 递归过程中存路
  std::vector<std::vector<int>> m_Answer[5];  // 长度3-7的环

  const int LIMIT_STEP = 3;
  std::unordered_set<int> m_StepArrive[MAXN];  // 一个点STEP步可以到达的结点

 private:
  struct Edge {
    int v, w;
  };
  std::vector<Edge> m_Sons[MAXN];               // 边
  int m_edgeNum = 0;                            // 边数目
  int m_dfn[MAXN], m_low[MAXN], m_stack[MAXN];  // Tarjan
  int m_catrgory[MAXN];                         // 所在联通分量id
  int m_inDegree[MAXN], m_outDegree[MAXN];      // 出入度
  bool m_inStack[MAXN];                         // Tarjan标记在不在栈内
  bool m_vis[MAXN];                             // 找环标记
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
  void PreRunPoint();

 private:
  inline void addEdge(int u, int v, int w);
  void tarjan(int u);
  void handelCircle(const int &dep);
  void findCircle(const int &ctg, int u, int dep);
  void sortAnswer();
  inline int GetID(int x);
  int GetMapID(int x);
  void pretreatTravel(int &st, int u, int dep);
};

void XJBG::Init() { m_TempPath = std::vector<int>(7, -1); }
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

void XJBG::pretreatTravel(int &st, int u, int dep) {
  if (dep > LIMIT_STEP) return;
  m_StepArrive[st].insert(u);
  for (auto &it : m_Sons[u]) {
    int v = it.v;
    if (m_vis[v] || m_catrgory[v] != m_catrgory[st]) continue;
    pretreatTravel(st, v, dep + 1);
  }
}
void XJBG::PreRunPoint() {
  ScopeTime t;

  for (auto &v : m_Circles) {
    pretreatTravel(v, v, 0);
  }

#ifdef LOCAL
  std::cerr << "@ PreRun #";
  t.LogTime();
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
      m_catrgory[cur] = u;
      tmp.emplace_back(cur);
    } while (cur != u);
    if (tmp.size() >= 3) {
      ++m_useScc;
      m_Circles.insert(m_Circles.end(), tmp.begin(), tmp.end());
    }
    ++m_scc;
  }
}

void XJBG::handelCircle(const int &dep) {
  std::vector<int> tmp;
  for (int i = 0; i < dep; ++i) {
    m_inCircle[m_TempPath[i]] = true;
    tmp.emplace_back(this->GetID(m_TempPath[i]));
  }
  m_Answer[dep - 3].emplace_back(tmp);
  ++m_answers;
  // if (m_answers % 10000 == 0) std::cerr << m_answers << "\n";
}

void XJBG::findCircle(const int &ctg, int u, int dep) {
  if (dep > 7 || m_vis[u] || m_catrgory[u] != ctg) return;
  if (dep > (7 - LIMIT_STEP) &&
      m_StepArrive[u].find(m_TempPath[0]) == m_StepArrive[u].end()) {
    return;
  }

  m_vis[u] = true;

  for (auto &it : m_Sons[u]) {
    int v = it.v;
    if (v == m_TempPath[0] && dep >= 3) {
      this->handelCircle(dep);
      continue;
    }
    m_TempPath[dep] = v;
    this->findCircle(ctg, v, dep + 1);
  }
  m_vis[u] = false;
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

void XJBG::FindPath() {
  ScopeTime t;

  for (auto &v : m_Circles) {
    if (m_inCircle[v] && m_inDegree[v] == m_outDegree[v] &&
        m_inDegree[v] == 1) {
      continue;
    }
    m_TempPath[0] = v;
    this->findCircle(m_catrgory[v], v, 1);
    m_vis[v] = true;
  }

#ifdef LOCAL
  std::cerr << "@ Answer: " << m_answers << "(";
  for (int i = 3; i < 7; ++i) {
    std::cerr << m_Answer[i - 3].size() << ",";
  }
  std::cerr << m_Answer[7 - 3].size() << ") #";
  t.LogTime();
#endif
}

void XJBG::sortAnswer() {
  for (auto &it : m_Answer) {
    std::sort(it.begin(), it.end());
  }

#ifdef TEST
  std::set<std::vector<int>> st;
  for (int i = 0; i < 5; ++i) {
    for (auto &it : m_Answer[i]) {
      std::vector<int> tmp;
      for (int k = 0; k < i + 3; ++k) {
        tmp.emplace_back(it[k]);
      }
      st.insert(tmp);
    }
  }
  std::cerr << "@ unique answer: " << st.size() << "\n";
#endif
}

void XJBG::SaveAnswer() {
  this->sortAnswer();
  std::ofstream fout(RESULT);
  fout << m_answers << "\n";
  for (int i = 0; i < 5; ++i) {
    for (auto &it : m_Answer[i]) {
      for (int k = 0; k < i + 3; ++k) {
        if (k != 0) fout << ",";
        fout << it[k];
      }
      fout << "\n";
    }
  }
  fout.close();
}

int main() {
  std::cerr << std::fixed << std::setprecision(3);

  XJBG *xjbg = new XJBG();
  xjbg->Init();
  xjbg->LoadData();
  xjbg->TarJan();
  xjbg->PreRunPoint();
  xjbg->FindPath();
  xjbg->SaveAnswer();
  return 0;
}