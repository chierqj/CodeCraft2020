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
#define TRAIN "../data/gen7w/test_data.txt"
#define RESULT "../data/gen7w/result.txt"
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
  static const int MAXN = 30000 + 7;
  static const int MAXM = 280000 + 7;

 private:
  struct Edge {
    int v, w;
  };
  std::vector<Edge> m_Edges[MAXN];
  int m_edgeNum = 0;
  int m_dfn[MAXN];
  int m_low[MAXN];
  int m_catrgory[MAXN];
  int m_stack[MAXN];
  bool m_inStack[MAXN], m_vis[MAXN], m_NotCircle[MAXN];
  int m_tarjanCount = 0, m_stackTop = 0, m_scc = 0, m_useScc = 0;

  std::vector<int> m_Circles;
  std::vector<int> m_TempPath;
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
  void findCircle(const int &ctg, int u, int dep);
  void sortAnswer();
};

void XJBG::Init() { m_TempPath = std::vector<int>(7, -1); }

inline void XJBG::addEdge(int u, int v, int w) {
  m_Edges[u].emplace_back(Edge{v, w});
  ++m_edgeNum;
}

void XJBG::LoadData() {
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
  std::cerr << "@ edge: " << m_edgeNum << "\n";
#endif
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

  for (auto &it : m_Edges[u]) {
    int v = it.v;
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

void XJBG::findCircle(const int &ctg, int u, int dep) {
  if (dep > 7 || m_vis[u] || m_catrgory[u] != ctg) return;
  m_vis[u] = true;

  for (auto &it : m_Edges[u]) {
    int v = it.v;

    if (v == m_TempPath[0] && dep >= 3) {
      m_Answer[dep - 3].emplace_back(m_TempPath);
      ++m_answers;
      continue;
    }

    m_TempPath[dep] = v;
    this->findCircle(ctg, v, dep + 1);
  }

  m_vis[u] = false;
}

void XJBG::TarJan() {
  for (int i = 0; i < MAXN; ++i) {
    if (!m_dfn[i] && !m_Edges[i].empty()) {
      this->tarjan(-1, i);
    }
  }
  std::sort(m_Circles.begin(), m_Circles.end());

#ifdef LOCAL
  std::cerr << "@ scc: " << m_scc << ", usescc: " << m_useScc << "\n";
#endif
}

void XJBG::FindPath() {
  ScopeTime t;
  for (auto &v : m_Circles) {
    m_TempPath[0] = v;
    this->findCircle(m_catrgory[v], v, 1);
    m_vis[v] = true;
  }

#ifdef LOCAL
  std::cerr << "@ total: " << m_answers << " (";
  for (int i = 3; i < 7; ++i) {
    std::cerr << m_Answer[i - 3].size() << ", ";
  }
  std::cerr << m_Answer[7 - 3].size() << ") ";
  t.LogTime();
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
  xjbg->FindPath();
  xjbg->SaveAnswer();
  return 0;
}