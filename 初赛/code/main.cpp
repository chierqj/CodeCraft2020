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
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#ifdef LOCAL
#define TRAIN "../data/std/test_data.txt"
#define RESULT "../data/std/result.txt"
#else
#define TRAIN "/data/test_data.txt"
#define RESULT "/projects/student/result.txt"
#endif
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
  static const int MAXN = 30000 + 7;
  static const int MAXM = 280000 + 7;

 private:
  struct Edge {
    int u, v, w;
    int next;
  } m_Edges[MAXM];
  int m_edgeNum = 0;
  int m_dfn[MAXN];
  int m_low[MAXN];
  int m_index[MAXN];
  int m_catrgory[MAXN];
  int m_stack[MAXN];
  bool m_inStack[MAXN], m_vis[MAXN];
  int m_tarjanCount = 0, m_stackTop = 0, m_scc = 0, m_useScc = 0;

  std::vector<int> m_Circles;
  std::vector<std::vector<int>> m_Answer[5];
  int m_answers = 0;

 public:
  void Init();
  void LoadData();
  void TarJan();
  void FindPath();
  void SaveAnswer();

 private:
  inline void addEdge(int u, int v, int w);
  void tarjan(int fa, int u);
  void findCircle(const int st, const int ctg, int u, int fa,
                  std::vector<int> path);
  void sortAnswer();
};

void XJBG::Init() { memset(m_index, -1, sizeof(m_index)); }

inline void XJBG::addEdge(int u, int v, int w) {
  m_Edges[m_edgeNum] = Edge{u, v, w, m_index[u]};
  m_index[u] = m_edgeNum++;
}

void XJBG::LoadData() {
  std::ifstream fin(TRAIN);
  int u, v, w;
  char c;
  while (fin >> u >> c >> v >> c >> w) {
    addEdge(u, v, w);
  }
  std::cerr << "@ edge: " << m_edgeNum << "\n";
  fin.close();

#ifdef TEST
  std::set<std::vector<int>> st;
  for (int i = 0; i < m_edgeNum; ++i) {
    st.insert({m_Edges[i].u, m_Edges[i].v});
  }
  std::cerr << "@ unique edge: " << st.size() << "\n";
#endif
}

void XJBG::tarjan(int fa, int u) {
  m_dfn[u] = m_low[u] = ++m_tarjanCount;
  m_stack[m_stackTop++] = u;
  m_inStack[u] = 1;
  for (int i = m_index[u]; i != -1; i = m_Edges[i].next) {
    int v = m_Edges[i].v;
    if (!m_dfn[v]) {
      this->tarjan(u, v);
      m_low[u] = std::min(m_low[u], m_low[v]);
    } else if (v != fa && m_inStack[v]) {
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

void XJBG::findCircle(const int st, const int ctg, int fa, int u,
                      std::vector<int> path) {
  int sz = path.size();
  if (sz > 7) return;

  for (int i = m_index[u]; i != -1; i = m_Edges[i].next) {
    int v = m_Edges[i].v;
    if (v == st && sz >= 3 && sz <= 7) {
      m_Answer[sz - 3].emplace_back(path);
      ++m_answers;
      continue;
    }

    if (v == fa || m_catrgory[v] != ctg || m_vis[v]) continue;

    path.emplace_back(v);
    m_vis[v] = true;
    this->findCircle(st, ctg, u, v, path);
    path.pop_back();
    m_vis[v] = false;
  }
}

void XJBG::TarJan() {
  for (int i = 0; i < MAXN; ++i) {
    if (!m_dfn[i] && m_index[i] != -1) {
      this->tarjan(-1, i);
    }
  }
  std::sort(m_Circles.begin(), m_Circles.end());

#ifdef LOCAL
  std::cerr << "@ scc: " << m_scc << ", usescc: " << m_useScc << "\n";
#endif
}

void XJBG::FindPath() {
  for (auto &v : m_Circles) {
    m_vis[v] = true;
    std::vector<int> path{v};
    this->findCircle(v, m_catrgory[v], -1, v, path);
  }

#ifdef LOCAL
  std::cerr << "@ total: " << m_answers << " (";
  for (int i = 3; i < 7; ++i) {
    std::cerr << m_Answer[i - 3].size() << ", ";
  }
  std::cerr << m_Answer[7 - 3].size() << ")\n";
#endif
}

void XJBG::sortAnswer() {
  for (auto &it : m_Answer) {
    std::sort(it.begin(), it.end());
  }

#ifdef TEST
  std::set<std::vector<int>> st;
  for (auto &ans : m_Answer) {
    for (auto &it : ans) {
      st.insert(it);
    }
  }
  std::cerr << "@ unique answer: " << st.size() << "\n";
#endif
}

void XJBG::SaveAnswer() {
  this->sortAnswer();
  std::ofstream fout(RESULT);
  fout << m_answers << "\n";
  for (auto &ans : m_Answer) {
    for (auto &it : ans) {
      for (int i = 0; i < it.size(); ++i) {
        if (i != 0) fout << ",";
        fout << it[i];
      }
      fout << "\n";
    }
  }
  fout.close();
  usleep(600000);
}

int main() {
  XJBG *xjbg = new XJBG();
  xjbg->Init();
  xjbg->LoadData();
  xjbg->TarJan();
  xjbg->FindPath();
  xjbg->SaveAnswer();
  return 0;
}
