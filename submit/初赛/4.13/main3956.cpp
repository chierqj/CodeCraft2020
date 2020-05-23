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
  inline void addEdge(int u, int v, int w);  // 加边
  void BackSearch(int &st, int dep);         // 反向剪枝
  void doFindCircle(int u, int dep);         // 正向找环
  void createAnswerBuffers(int pid);         // 解析int
  inline void connectBuffer(char *buffer, int &node, int &idx, char c);

 public:
  static const int MAXN = 200000 + 7;                   // 总点数
  static const int NTHREAD = 8;                         // 线程个数
  const int PARAM[NTHREAD] = {1, 2, 3, 4, 5, 6, 7, 8};  // 线程权重
  // const int PARAM[NTHREAD] = {2, 3, 5, 40};  // 线程权重

 private:
  int m_maxID = 0;                              // 最大点
  int m_Edges[MAXN][3][100];                    // u->[son, father, weight]
  int m_CountSons[MAXN], m_CountFathers[MAXN];  // 边数目
  int m_edgeNum = 0;                            // 边数目
  std::vector<PreBuffer> m_MapID;               // 预处理答案
  std::vector<int> m_Circles;                   // 点集合

 private:
  int m_answers = 0;                   // 环个数
  int m_start = 0;                     // 起点
  std::vector<bool> m_OneReachable;    // 一步可达
  std::vector<bool> m_TwoReachable;    // 两步可达
  std::vector<bool> m_ThreeReachable;  // 三步可达
  Vector<char *> m_Answer[5];          // 结果
  Vector<int> m_Count[5];              // 答案个数
  char *m_AnsBuf[5];                   // 答案
  uint32_t m_AnsBufLen[5];             // bufferLength

 private:
  int m_TotalAnswers = shmget(IPC_PRIVATE, sizeof(int), IPC_CREAT | 0600);
  int m_FlagFork = shmget(IPC_PRIVATE, NTHREAD * sizeof(int), IPC_CREAT | 0600);
  int m_FlagSave = shmget(IPC_PRIVATE, sizeof(int), IPC_CREAT | 0600);
  uint32_t m_TotalBufferSize =
      shmget(IPC_PRIVATE, sizeof(uint32_t), IPC_CREAT | 0600);
  uint32_t m_BufferStart =
      shmget(IPC_PRIVATE, sizeof(uint32_t), IPC_CREAT | 0600);

  int *m_TotalAnswersPtr;
  int *m_FlagForkPtr;
  int *m_FlagSavePtr;
  uint32_t *m_TotalBufferSizePtr;
  uint32_t *m_BufferStartPtr;
};

void XJBG::GetShmPtr() {
  m_TotalAnswersPtr = (int *)shmat(m_TotalAnswers, NULL, 0);
  m_FlagForkPtr = (int *)shmat(m_FlagFork, NULL, 0);
  m_FlagSavePtr = (int *)shmat(m_FlagSave, NULL, 0);
  m_TotalBufferSizePtr = (uint32_t *)shmat(m_TotalBufferSize, NULL, 0);
  m_BufferStartPtr = (uint32_t *)shmat(m_BufferStart, NULL, 0);
}

inline void XJBG::addEdge(int u, int v, int w) {
  m_Edges[u][0][m_CountSons[u]] = v;
  m_Edges[u][2][m_CountSons[u]++] = w;
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

void XJBG::BackSearch(int &st, int dep) {
  for (int i = 0; i < m_CountFathers[st]; ++i) {
    int v1 = m_Edges[st][1][i];
    if (v1 <= st) continue;
    m_OneReachable[v1] = true;
    m_TwoReachable[v1] = true;
    m_ThreeReachable[v1] = true;
    for (int j = 0; j < m_CountFathers[v1]; ++j) {
      int v2 = m_Edges[v1][1][j];
      if (v2 <= st || v2 == v1) continue;
      m_TwoReachable[v2] = true;
      m_ThreeReachable[v2] = true;
      for (int k = 0; k < m_CountFathers[v2]; ++k) {
        int v3 = m_Edges[v2][1][k];
        if (v3 <= st || v3 == v2 || v3 == v1) continue;
        m_ThreeReachable[v3] = true;
      }
    }
  }
}

inline void XJBG::connectBuffer(char *buffer, int &node, int &idx, char c) {
  memcpy(buffer + idx, m_MapID[node].str + m_MapID[node].start,
         m_MapID[node].length);
  idx += m_MapID[node].length;
  buffer[idx++] = c;
}
void XJBG::doFindCircle(int u, int dep) {
  int sz0 = m_MapID[u].length;
  for (int it1 = 0; it1 < m_CountSons[u]; ++it1) {
    int v1 = m_Edges[u][0][it1];
    if (v1 <= u) continue;
    for (int it2 = 0; it2 < m_CountSons[v1]; ++it2) {
      int v2 = m_Edges[v1][0][it2];
      if (v2 <= u || v2 == v1) continue;
      for (int it3 = 0; it3 < m_CountSons[v2]; ++it3) {
        int v3 = m_Edges[v2][0][it3];
        if (v3 == u) {
          int sz = m_MapID[u].length + m_MapID[v1].length + m_MapID[v2].length;
          char *buf = new char[sz + 3];
          int idx = 0;
          connectBuffer(buf, u, idx, ',');
          connectBuffer(buf, v1, idx, ',');
          connectBuffer(buf, v2, idx, '\n');
          m_Answer[0].emplace_back(buf);
          m_Count[0].emplace_back(idx);
          ++m_answers;
          continue;
        }
        if (v3 <= u || v3 == v2 || v3 == v1) continue;
        for (int it4 = 0; it4 < m_CountSons[v3]; ++it4) {
          int v4 = m_Edges[v3][0][it4];
          if (v4 == u) {
            int sz = m_MapID[u].length + m_MapID[v1].length +
                     m_MapID[v2].length + m_MapID[v3].length;

            char *buf = new char[sz + 4];
            int idx = 0;
            connectBuffer(buf, u, idx, ',');
            connectBuffer(buf, v1, idx, ',');
            connectBuffer(buf, v2, idx, ',');
            connectBuffer(buf, v3, idx, '\n');
            m_Answer[1].emplace_back(buf);
            m_Count[1].emplace_back(idx);
            ++m_answers;
            continue;
          }
          if (v4 <= u || v4 == v3 || v4 == v2 || v4 == v1 ||
              !m_ThreeReachable[v4])
            continue;
          for (int it5 = 0; it5 < m_CountSons[v4]; ++it5) {
            int v5 = m_Edges[v4][0][it5];
            if (v5 == u) {
              int sz = m_MapID[u].length + m_MapID[v1].length +
                       m_MapID[v2].length + m_MapID[v3].length +
                       m_MapID[v4].length;
              char *buf = new char[sz + 5];
              int idx = 0;
              connectBuffer(buf, u, idx, ',');
              connectBuffer(buf, v1, idx, ',');
              connectBuffer(buf, v2, idx, ',');
              connectBuffer(buf, v3, idx, ',');
              connectBuffer(buf, v4, idx, '\n');
              m_Answer[2].emplace_back(buf);
              m_Count[2].emplace_back(idx);
              ++m_answers;
              continue;
            }
            if (v5 <= u || v5 == v4 || v5 == v3 || v5 == v2 || v5 == v1 ||
                !m_TwoReachable[v5])
              continue;
            for (int it6 = 0; it6 < m_CountSons[v5]; ++it6) {
              int v6 = m_Edges[v5][0][it6];
              if (v6 == u) {
                int sz = m_MapID[u].length + m_MapID[v1].length +
                         m_MapID[v2].length + m_MapID[v3].length +
                         m_MapID[v4].length + m_MapID[v5].length;
                char *buf = new char[sz + 6];
                int idx = 0;
                connectBuffer(buf, u, idx, ',');
                connectBuffer(buf, v1, idx, ',');
                connectBuffer(buf, v2, idx, ',');
                connectBuffer(buf, v3, idx, ',');
                connectBuffer(buf, v4, idx, ',');
                connectBuffer(buf, v5, idx, '\n');
                m_Answer[3].emplace_back(buf);
                m_Count[3].emplace_back(idx);
                ++m_answers;
                continue;
              }
              if (v6 <= u || v6 == v5 || v6 == v4 || v6 == v3 || v6 == v2 ||
                  v6 == v1)
                continue;
              if (m_OneReachable[v6]) {
                int sz = m_MapID[u].length + m_MapID[v1].length +
                         m_MapID[v2].length + m_MapID[v3].length +
                         m_MapID[v4].length + m_MapID[v5].length +
                         m_MapID[v6].length;
                char *buf = new char[sz + 7];
                int idx = 0;
                connectBuffer(buf, u, idx, ',');
                connectBuffer(buf, v1, idx, ',');
                connectBuffer(buf, v2, idx, ',');
                connectBuffer(buf, v3, idx, ',');
                connectBuffer(buf, v4, idx, ',');
                connectBuffer(buf, v5, idx, ',');
                connectBuffer(buf, v6, idx, '\n');
                m_Answer[4].emplace_back(buf);
                m_Count[4].emplace_back(idx);
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
      std::swap(m_Edges[v][2][j], m_Edges[v][2][minIndex]);
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

  for (int i = st; i < ed; ++i) {
    m_OneReachable = std::vector<bool>(m_maxID, false);
    m_TwoReachable = std::vector<bool>(m_maxID, false);
    m_ThreeReachable = std::vector<bool>(m_maxID, false);
    int v = m_Circles[i];
    m_start = v;
    this->BackSearch(v, 0);
    this->doFindCircle(v, 1);
  }
  *m_TotalAnswersPtr += m_answers;
  this->createAnswerBuffers(pid);
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
void XJBG::createAnswerBuffers(int pid) {
  for (int len = 0; len < 5; ++len) {
    auto &ans = m_Answer[len];
    uint32_t tbufsize = (uint32_t)ans.n * (len + 3) * 11;
    m_AnsBuf[len] = new char[tbufsize];
    uint32_t tidx = 0;
    for (int k = 0; k < ans.n; ++k) {
      memcpy(m_AnsBuf[len] + tidx, ans[k], m_Count[len][k]);
      tidx += m_Count[len][k];
    }
    m_AnsBufLen[len] = tidx;
    *m_TotalBufferSizePtr += tidx;
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

    int fd = open(RESULT, O_RDWR | O_CREAT, 0666);
    uint32_t bufsize = *m_TotalBufferSizePtr + (12 - firIdx);
    ftruncate(fd, bufsize);
    char *result =
        (char *)mmap(NULL, bufsize, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    close(fd);
    memcpy(result, firBuf + firIdx, 12 - firIdx);
    *m_BufferStartPtr = 12 - firIdx;
    *m_FlagSavePtr = 0;
  }

  for (int len = 0; len < 5; ++len) {
    while (*m_FlagSavePtr != pid * 5 + len) usleep(1);
    struct stat sb;
    int fd = open(RESULT, O_RDWR | O_CREAT, 0666);
    fstat(fd, &sb);
    char *result = (char *)mmap(NULL, sb.st_size, PROT_READ | PROT_WRITE,
                                MAP_SHARED, fd, 0);
    close(fd);
    memcpy(result + *m_BufferStartPtr, m_AnsBuf[len], m_AnsBufLen[len]);
    *m_BufferStartPtr += m_AnsBufLen[len];
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
    std::cerr << "TotalBufferSize: " << *m_TotalBufferSizePtr << "\n";
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

  // for (int i = 0; i < XJBG::NTHREAD - 1; ++i) {
  //   waitpid(Children[i], NULL, 0);
  // }
  return 0;
}