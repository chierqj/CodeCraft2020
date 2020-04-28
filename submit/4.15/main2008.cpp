#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/ipc.h>
#include <sys/mman.h>
#include <sys/shm.h>
#include <sys/stat.h>
#include <sys/types.h>
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
#include <random>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <thread>
#include <tuple>
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
template <typename T>
struct Vector {
  T *a;
  int n;
  Vector() {
    a = nullptr;
    n = 0;
  }
  void emplace_back(T x) {
    if ((n & -n) == n) {
      a = (T *)realloc(a, (n * 2 + 1) * sizeof(T));
    }
    a[n++] = x;
  }
  void init(int num, T val) {
    a = nullptr;
    n = 0;
    for (int i = 0; i < num; ++i) {
      if ((n & -n) == n) {
        a = (T *)realloc(a, (n * 2 + 1) * sizeof(T));
      }
      a[n++] = val;
    }
  }
  T &operator[](int idx) { return a[idx]; }
  // struct Iterator {
  //   int index;
  //   Vector &outer;
  //   Iterator(Vector &o, int i) : outer(o), index(i) {}
  //   void operator++() { index++; }
  //   T operator*() const { return outer[index]; }
  //   bool operator!=(const Iterator &i) { return i.index != index; }
  // };
  // Iterator begin() { return Iterator(*this, 0); }
  // Iterator end() { return Iterator(*this, n); }
};

class XJBG {
  struct PreBuffer {
    char str[10];
    int start;
    int length;
  };

 private:
  std::string m_debug[8] = {
      "",           "\t\t",         "\t\t\t",         "\t\t\t\t",
      "\t\t\t\t\t", "\t\t\t\t\t\t", "\t\t\t\t\t\t\t", "\t\t\t\t\t\t\t\t"};

 public:
  void LoadData();           // 加载数据
  void FindCircle(int pid);  // 找环
  void PreSave();            // 预处理ID
  void SaveAnswer(int pid);  // 保存答案
  void SortEdge();           // 儿子排序
  void WaitFork();           // 等待进程
  void GetShmPtr();          // 获取共享数据句柄

 private:
  inline void addEdge(int u, int v, int w);              // 加边
  void BackSearch(int st, int dep, Vector<int> &tmp);    // 反向剪枝
  void doFindCircle(int st);                             // 正向找环
  inline void connectBuffer(int dep, int node, char c);  // dfs拼buffer

 public:
  static const int MAXN = 200000 + 7;  // 总点数
  static const int NTHREAD = 16;       // 线程个数
  const int PARAM[NTHREAD] = {1, 2,  3,  4,  5,  6,  7,  8,
                              9, 10, 11, 12, 13, 14, 15, 16};  // 线程权重
  // const int PARAM[NTHREAD] = {2, 3, 5, 40};  // 线程权重

 private:
  int m_maxID = 0;                              // 最大点
  int m_Edges[MAXN][2][50];                     // u->[son, father, weight]
  int m_CountSons[MAXN], m_CountFathers[MAXN];  // 边数目
  int m_edgeNum = 0;                            // 边数目
  std::vector<PreBuffer> m_MapID;               // 预处理答案
  std::vector<int> m_Circles;                   // 点集合
  int m_answers = 0;                            // 环个数
  Vector<char> m_Reachable;                     // 可达
  char *m_Answer[5];                            // 答案
  uint32_t m_AnswerLength[5];                   // bufferLength

 private:
  int m_TotalAnswers = shmget(IPC_PRIVATE, sizeof(int), IPC_CREAT | 0600);
  int m_FlagFork = shmget(IPC_PRIVATE, NTHREAD * sizeof(int), IPC_CREAT | 0600);
  int m_FlagSave = shmget(IPC_PRIVATE, sizeof(int), IPC_CREAT | 0600);

  int *m_TotalAnswersPtr;
  int *m_FlagForkPtr;
  int *m_FlagSavePtr;
};

void XJBG::GetShmPtr() {
  m_TotalAnswersPtr = (int *)shmat(m_TotalAnswers, NULL, 0);
  m_FlagForkPtr = (int *)shmat(m_FlagFork, NULL, 0);
  m_FlagSavePtr = (int *)shmat(m_FlagSave, NULL, 0);
}

inline void XJBG::addEdge(int u, int v, int w) {
  m_Edges[u][0][m_CountSons[u]++] = v;
  // m_Edges[u][2][m_CountSons[u]++] = w;
  m_Edges[v][1][m_CountFathers[v]++] = u;
  ++m_edgeNum;
  m_maxID = std::max(m_maxID, std::max(u, v) + 1);
}

void XJBG::LoadData() {
  ScopeTime t;

  struct stat sb;
  int fd = open(TRAIN, O_RDONLY);
  fstat(fd, &sb);
  long long bufsize = sb.st_size;
  char *buffer = (char *)mmap(NULL, bufsize, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);

  int u = 0, v = 0, w = 0;
  char *ptr = buffer;
  while (ptr - buffer < bufsize) {
    while (*ptr != ',') {
      u = u * 10 + *ptr - '0';
      ++ptr;
    }
    ++ptr;
    while (*ptr != ',') {
      v = v * 10 + *ptr - '0';
      ++ptr;
    }
    ++ptr;
    while (*ptr != '\n') {
      w = w * 10 + *ptr - '0';
      ++ptr;
    }
    ++ptr;
    addEdge(u, v, w);
    u = v = w = 0;
  }

  for (int i = 0; i < m_maxID; ++i) {
    if (m_CountSons[i] > 0 && m_CountFathers[i] > 0) {
      m_Circles.emplace_back(i);
    }
  }

#ifdef TEST
  std::cerr << "@ LoadData: (V: " << m_maxID << ", E: " << m_edgeNum << ") #";
  t.LogTime();
#endif
}

void XJBG::BackSearch(int st, int dep, Vector<int> &tmp) {
  for (int i = 0; i < m_CountFathers[st]; ++i) {
    int v1 = m_Edges[st][1][i];
    if (v1 <= st) continue;
    m_Reachable[v1] |= 7;
    tmp.emplace_back(v1);
    for (int j = 0; j < m_CountFathers[v1]; ++j) {
      int v2 = m_Edges[v1][1][j];
      if (v2 <= st || v2 == v1) continue;
      m_Reachable[v2] |= 6;
      tmp.emplace_back(v2);
      for (int k = 0; k < m_CountFathers[v2]; ++k) {
        int v3 = m_Edges[v2][1][k];
        if (v3 <= st || v3 == v2 || v3 == v1) continue;
        m_Reachable[v3] |= 4;
        tmp.emplace_back(v3);
      }
    }
  }
}

inline void XJBG::connectBuffer(int dep, int node, char c) {
  const auto &mpid = m_MapID[node];
  memcpy(m_Answer[dep] + m_AnswerLength[dep], mpid.str + mpid.start,
         mpid.length);
  m_AnswerLength[dep] += mpid.length;
  m_Answer[dep][m_AnswerLength[dep]++] = c;
}

void XJBG::doFindCircle(int st) {
  for (int it1 = 0; it1 < m_CountSons[st]; ++it1) {
    int v1 = m_Edges[st][0][it1];
    if (v1 <= st) continue;
    for (int it2 = 0; it2 < m_CountSons[v1]; ++it2) {
      int v2 = m_Edges[v1][0][it2];
      if (v2 <= st || v2 == v1) continue;
      for (int it3 = 0; it3 < m_CountSons[v2]; ++it3) {
        int v3 = m_Edges[v2][0][it3];
        if (v3 == st) {
          connectBuffer(0, st, ',');
          connectBuffer(0, v1, ',');
          connectBuffer(0, v2, '\n');
          ++m_answers;
          continue;
        }
        if (v3 <= st || v3 == v2 || v3 == v1) continue;
        for (int it4 = 0; it4 < m_CountSons[v3]; ++it4) {
          int v4 = m_Edges[v3][0][it4];
          if (v4 == st) {
            connectBuffer(1, st, ',');
            connectBuffer(1, v1, ',');
            connectBuffer(1, v2, ',');
            connectBuffer(1, v3, '\n');
            ++m_answers;
            continue;
          }
          if (v4 <= st || v4 == v3 || v4 == v2 || v4 == v1 ||
              !(m_Reachable[v4] & 4))
            continue;
          for (int it5 = 0; it5 < m_CountSons[v4]; ++it5) {
            int v5 = m_Edges[v4][0][it5];
            if (v5 == st) {
              connectBuffer(2, st, ',');
              connectBuffer(2, v1, ',');
              connectBuffer(2, v2, ',');
              connectBuffer(2, v3, ',');
              connectBuffer(2, v4, '\n');
              ++m_answers;
              continue;
            }
            if (v5 <= st || v5 == v4 || v5 == v3 || v5 == v2 || v5 == v1 ||
                !(m_Reachable[v5] & 2))
              continue;
            for (int it6 = 0; it6 < m_CountSons[v5]; ++it6) {
              int v6 = m_Edges[v5][0][it6];
              if (v6 == st) {
                connectBuffer(3, st, ',');
                connectBuffer(3, v1, ',');
                connectBuffer(3, v2, ',');
                connectBuffer(3, v3, ',');
                connectBuffer(3, v4, ',');
                connectBuffer(3, v5, '\n');
                ++m_answers;
                continue;
              }
              if (v6 <= st || v6 == v5 || v6 == v4 || v6 == v3 || v6 == v2 ||
                  v6 == v1)
                continue;
              if (m_Reachable[v6] & 1) {
                connectBuffer(4, st, ',');
                connectBuffer(4, v1, ',');
                connectBuffer(4, v2, ',');
                connectBuffer(4, v3, ',');
                connectBuffer(4, v4, ',');
                connectBuffer(4, v5, ',');
                connectBuffer(4, v6, '\n');
                ++m_answers;
              }
            }
          }
        }
      }
    }
  }
}

void XJBG::SortEdge() {
  for (auto &v : m_Circles) {
    int count = m_CountSons[v];
    for (int j = 0; j < count; ++j) {
      int minIndex = j;
      for (int k = minIndex + 1; k < count; ++k) {
        if (m_Edges[v][0][k] < m_Edges[v][0][minIndex]) {
          minIndex = k;
        }
      }
      std::swap(m_Edges[v][0][j], m_Edges[v][0][minIndex]);
      // std::swap(m_Edges[v][2][j], m_Edges[v][2][minIndex]);
    }
  }
}

void XJBG::FindCircle(int pid) {
  ScopeTime t;

#ifdef TEST
  std::cerr << m_debug[pid] << pid << ": FindStart\n";
#endif
  int sum = 0;
  for (int i = 0; i < NTHREAD; ++i) sum += PARAM[i];
  int block = m_Circles.size() / sum, st = 0;
  for (int i = 0; i < pid; ++i) st += PARAM[i] * block;
  int ed = (pid == NTHREAD - 1 ? m_Circles.size() : st + block * PARAM[pid]);
  m_Reachable.init(m_maxID, 0);

  for (int i = 0; i < 5; ++i) {
    m_Answer[i] = new char[3000000 * (i + 3) * 11];
    m_AnswerLength[i] = 0;
  }

  for (int i = st; i < ed; ++i) {
    Vector<int> tmp;
    int v = m_Circles[i];
    this->BackSearch(v, 0, tmp);
    this->doFindCircle(v);
    for (int j = 0; j < tmp.n; ++j) {
      m_Reachable[tmp[j]] = 0;
    }
  }

  *m_TotalAnswersPtr += m_answers;
  m_FlagForkPtr[pid] = -1;
#ifdef TEST
  std::cerr << m_debug[pid] << pid << ": FindEnd\n";
  std::cerr << m_debug[pid] << pid << ": " << m_answers << "\n";
  std::cerr << m_debug[pid] << pid << ": ";
  t.LogTime();
#endif
}

void XJBG::PreSave() {
  m_MapID.resize(m_maxID);
  for (auto &v : m_Circles) {
    int x = v;
    if (x == 0) {
      m_MapID[v].start = 9;
      m_MapID[v].str[9] = '0';
      m_MapID[v].length = 1;
      continue;
    }
    int idx = 10;
    while (x) {
      m_MapID[v].str[--idx] = x % 10 + '0';
      x /= 10;
    }
    m_MapID[v].start = idx;
    m_MapID[v].length = 10 - idx;
  }
}
void XJBG::SaveAnswer(int pid) {
#ifdef TEST
  std::cerr << m_debug[pid] << pid << ": SaveStart\n";
#endif
  int answers = *m_TotalAnswersPtr;
  if (pid == 0) {
    int firIdx = 12;
    char firBuf[12];
    firBuf[--firIdx] = '\n';
    if (answers == 0) {
      firBuf[--firIdx] = '0';
    } else {
      while (answers) {
        firBuf[--firIdx] = answers % 10 + '0';
        answers /= 10;
      }
    }

    FILE *fp = fopen(RESULT, "w");
    fwrite(firBuf + firIdx, 1, 12 - firIdx, fp);
    fclose(fp);
    *m_FlagSavePtr = 0;
  }

  for (int len = 0; len < 5; ++len) {
    while (*m_FlagSavePtr != pid * 5 + len) usleep(1);

    FILE *fp = fopen(RESULT, "at");
    fwrite(m_Answer[len], 1, m_AnswerLength[len], fp);
    fclose(fp);

    if (pid == NTHREAD - 1) {
      *m_FlagSavePtr = len + 1;
    } else {
      *m_FlagSavePtr = (pid + 1) * 5 + len;
    }
  }

#ifdef TEST
  std::cerr << m_debug[pid] << pid << ": SaveEnd\n";
#endif
#ifdef LOCAL
  if (pid == 0) {
    std::cerr << "TotalAnswer: " << *m_TotalAnswersPtr << "\n";
  }
#endif
}

void XJBG::WaitFork() {
  while (true) {
    int x = 0;
    for (int i = 0; i < NTHREAD; ++i) x += m_FlagForkPtr[i];
    if (x == -NTHREAD) {
      return;
    }
    usleep(1);
  }
}

int main() {
  std::cerr << std::fixed << std::setprecision(3);

  XJBG *xjbg = new XJBG();
  xjbg->LoadData();
  xjbg->PreSave();
  xjbg->SortEdge();

  pid_t Children[XJBG::NTHREAD - 1] = {0};
  int pid = 0;
  for (int i = 0; i < XJBG::NTHREAD - 1; i++) {
    if (pid == 0) {
      int x = fork();
      Children[i] = x;
      if (x == 0) {
        pid = i + 1;
        break;
      }
    }
  }

  xjbg->GetShmPtr();
  xjbg->FindCircle(pid);
  xjbg->WaitFork();
  xjbg->SaveAnswer(pid);

  if (pid != 0) exit(0);
  return 0;
}