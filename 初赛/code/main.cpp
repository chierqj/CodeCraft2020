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
 public:
#ifdef XJBGJUDGE
  ~XJBG() {
    delete m_Answer0;
    delete m_Answer1;
    delete m_Answer2;
    delete m_Answer3;
    delete m_Answer4;
    delete m_AnswerLength0;
    delete m_AnswerLength1;
    delete m_AnswerLength2;
    delete m_AnswerLength3;
    delete m_AnswerLength4;
    delete m_ThreeCircle;
  }
#endif

 private:
  struct PreBuffer {
    char str[11];
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
  inline void addEdge(int u, int v, int w);            // 加边
  void BackSearch(int st, int dep, Vector<int> &tmp);  // 反向剪枝
  void doFindCircle(int st);                           // 正向找环
  inline void connectBuffer(char *answer, uint32_t *length, int node);

 public:
  static const int MAXN = 200000 + 7;  // 总点数
  static const int NTHREAD = 20;       // 线程个数
  const int PARAM[NTHREAD] = {1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11,
                              12, 13, 14, 15, 16, 17, 18, 19, 20};  // 线程权重
  const int PARAM_SUM = 210;

 private:
  int m_maxID = 0;  // 最大点
  // int m_Edges[2][MAXN][50];                     // u->[son, father, weight]
  int m_Sons[MAXN][50];                         //
  int m_Fathers[MAXN][50];                      //
  int m_CountSons[MAXN], m_CountFathers[MAXN];  // 边数目
  int m_edgeNum = 0;                            // 边数目
  std::vector<PreBuffer> m_MapID;               // 预处理答案
  Vector<int> m_Circles;                        // 点集合
  int m_answers = 0;                            // 环个数
  Vector<char> m_Reachable;                     // 可达
  char *m_Answer0;                              // 长度为3的环
  char *m_Answer1;                              // 长度为4的环
  char *m_Answer2;                              // 长度为5的环
  char *m_Answer3;                              // 长度为6的环
  char *m_Answer4;                              // 长度为7的环
  uint32_t *m_AnswerLength0;                    // 长度为3的环buffer长度
  uint32_t *m_AnswerLength1;                    // 长度为4的环buffer长度
  uint32_t *m_AnswerLength2;                    // 长度为5的环buffer长度
  uint32_t *m_AnswerLength3;                    // 长度为6的环buffer长度
  uint32_t *m_AnswerLength4;                    // 长度为7的环buffer长度
  char *m_ThreeCircle = new char[33];

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
  // m_Edges[0][u][m_CountSons[u]++] = v;
  // m_Edges[u][2][m_CountSons[u]++] = w;
  // m_Edges[1][v][m_CountFathers[v]++] = u;
  m_Sons[u][m_CountSons[u]++] = v;
  m_Fathers[v][m_CountFathers[v]++] = u;
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
  int v1, v2, v3;
  int count0 = m_CountFathers[st], count1, count2;
  m_Reachable[st] = 7;
  tmp.emplace_back(st);
  for (int i = 0; i < count0; ++i) {
    v1 = m_Fathers[st][i];
    if (v1 <= st) continue;
    m_Reachable[v1] |= 7;
    tmp.emplace_back(v1);
    count1 = m_CountFathers[v1];
    for (int j = 0; j < count1; ++j) {
      v2 = m_Fathers[v1][j];
      if (v2 <= st) continue;
      m_Reachable[v2] |= 6;
      tmp.emplace_back(v2);
      count2 = m_CountFathers[v2];
      for (int k = 0; k < count2; ++k) {
        v3 = m_Fathers[v2][k];
        if (v3 <= st || v3 == v1) continue;
        m_Reachable[v3] |= 4;
        tmp.emplace_back(v3);
      }
    }
  }
}

inline void XJBG::connectBuffer(char *answer, uint32_t *length, int node) {
  const auto &mpid = m_MapID[node];
  memcpy(answer + (*length), mpid.str + mpid.start, mpid.length);
  (*length) += mpid.length;
}

void XJBG::doFindCircle(int st) {
  const auto &mpid0 = m_MapID[st];
  memcpy(m_ThreeCircle, mpid0.str + mpid0.start, mpid0.length);
  int IDX0 = mpid0.length, IDX1 = 0, IDX2 = 0;

  int v1, v2, v3, v4, v5, v6, v7;
  int count0 = m_CountSons[st], count1, count2, count3, count4, count5;
  int *son0 = m_Sons[st], *son1, *son2, *son3, *son4, *son5;

  for (int it1 = 0; it1 < count0; ++it1) {
    v1 = *(son0 + it1);
    if (v1 < st) continue;

    const auto &mpid1 = m_MapID[v1];
    memcpy(m_ThreeCircle + IDX0, mpid1.str + mpid1.start, mpid1.length);
    IDX1 = IDX0 + mpid1.length;

    count1 = m_CountSons[v1];
    son1 = m_Sons[v1];
    for (int it2 = 0; it2 < count1; ++it2) {
      v2 = *(son1 + it2);
      if (v2 <= st) continue;

      const auto &mpid2 = m_MapID[v2];
      memcpy(m_ThreeCircle + IDX1, mpid2.str + mpid2.start, mpid2.length);
      IDX2 = IDX1 + mpid2.length;

      count2 = m_CountSons[v2];
      son2 = m_Sons[v2];
      for (int it3 = 0; it3 < count2; ++it3) {
        v3 = *(son2 + it3);
        if (v3 < st) {
          continue;
        } else if (!(v3 ^ st)) {
          memcpy(m_Answer0 + (*m_AnswerLength0), m_ThreeCircle, IDX2);
          (*m_AnswerLength0) += IDX2;
          *(m_Answer0 + (*m_AnswerLength0) - 1) = '\n';
          ++m_answers;
          continue;
        } else if (!(v3 ^ v1)) {
          continue;
        }
        count3 = m_CountSons[v3];
        son3 = m_Sons[v3];
        for (int it4 = 0; it4 < count3; ++it4) {
          v4 = *(son3 + it4);
          if (!(m_Reachable[v4] & 4)) {
            continue;
          } else if (!(v4 ^ st)) {
            memcpy(m_Answer1 + (*m_AnswerLength1), m_ThreeCircle, IDX2);
            (*m_AnswerLength1) += IDX2;
            connectBuffer(m_Answer1, m_AnswerLength1, v3);
            *(m_Answer1 + (*m_AnswerLength1) - 1) = '\n';
            ++m_answers;
            continue;
          } else if (!(v4 ^ v2) || !(v4 ^ v1)) {
            continue;
          }
          count4 = m_CountSons[v4];
          son4 = m_Sons[v4];
          for (int it5 = 0; it5 < count4; ++it5) {
            v5 = *(son4 + it5);
            if (!(m_Reachable[v5] & 2)) {
              continue;
            } else if (!(v5 ^ st)) {
              memcpy(m_Answer2 + (*m_AnswerLength2), m_ThreeCircle, IDX2);
              (*m_AnswerLength2) += IDX2;
              connectBuffer(m_Answer2, m_AnswerLength2, v3);
              connectBuffer(m_Answer2, m_AnswerLength2, v4);
              *(m_Answer2 + (*m_AnswerLength2) - 1) = '\n';
              ++m_answers;
              continue;
            } else if (!(v5 ^ v3) || !(v5 ^ v2) || !(v5 ^ v1)) {
              continue;
            }
            count5 = m_CountSons[v5];
            son5 = m_Sons[v5];
            for (int it6 = 0; it6 < count5; ++it6) {
              v6 = *(son5 + it6);
              if (!(m_Reachable[v6] & 1)) {
                continue;
              } else if (!(v6 ^ st)) {
                memcpy(m_Answer3 + (*m_AnswerLength3), m_ThreeCircle, IDX2);
                (*m_AnswerLength3) += IDX2;
                connectBuffer(m_Answer3, m_AnswerLength3, v3);
                connectBuffer(m_Answer3, m_AnswerLength3, v4);
                connectBuffer(m_Answer3, m_AnswerLength3, v5);
                *(m_Answer3 + (*m_AnswerLength3) - 1) = '\n';
                ++m_answers;
                continue;
              } else if (!(v6 ^ v4) || !(v6 ^ v3) || !(v6 ^ v2) || !(v6 ^ v1)) {
                continue;
              }
              memcpy(m_Answer4 + (*m_AnswerLength4), m_ThreeCircle, IDX2);
              (*m_AnswerLength4) += IDX2;
              connectBuffer(m_Answer4, m_AnswerLength4, v3);
              connectBuffer(m_Answer4, m_AnswerLength4, v4);
              connectBuffer(m_Answer4, m_AnswerLength4, v5);
              connectBuffer(m_Answer4, m_AnswerLength4, v6);
              *(m_Answer4 + (*m_AnswerLength4) - 1) = '\n';
              ++m_answers;
            }
          }
        }
      }
    }
  }
}

void XJBG::SortEdge() {
  for (int i = 0; i < m_Circles.n; ++i) {
    int v = m_Circles[i];
    int count = m_CountSons[v];
    for (int j = 0; j < count; ++j) {
      int minIndex = j;
      for (int k = minIndex + 1; k < count; ++k) {
        if (m_Sons[v][k] < m_Sons[v][minIndex]) {
          minIndex = k;
        }
      }
      std::swap(m_Sons[v][j], m_Sons[v][minIndex]);
    }
  }
}

void XJBG::FindCircle(int pid) {
  ScopeTime t;

#ifdef TEST
  std::cerr << m_debug[pid] << pid << ": FindStart\n";
#endif
  int block = m_Circles.n / PARAM_SUM, st = 0;
  for (int i = 0; i < pid; ++i) st += PARAM[i] * block;
  int ed = (pid == NTHREAD - 1 ? m_Circles.n : st + block * PARAM[pid]);
  m_Reachable.init(m_maxID, 0);

  m_Answer0 = new char[3000000 * 33];
  m_Answer1 = new char[3000000 * 44];
  m_Answer2 = new char[3000000 * 55];
  m_Answer3 = new char[3000000 * 66];
  m_Answer4 = new char[3000000 * 77];
  m_AnswerLength0 = new uint32_t;
  m_AnswerLength1 = new uint32_t;
  m_AnswerLength2 = new uint32_t;
  m_AnswerLength3 = new uint32_t;
  m_AnswerLength4 = new uint32_t;

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
  m_FlagForkPtr[pid] = 1;
#ifdef TEST
  std::cerr << m_debug[pid] << pid << ": FindEnd\n";
  std::cerr << m_debug[pid] << pid << ": " << m_answers << "\n";
  std::cerr << m_debug[pid] << pid << ": ";
  t.LogTime();
#endif
}

void XJBG::PreSave() {
  m_MapID.resize(m_maxID);
  for (int i = 0; i < m_Circles.n; ++i) {
    int v = m_Circles[i];
    int x = v;
    if (x == 0) {
      m_MapID[v].start = 9;
      m_MapID[v].str[10] = ',';
      m_MapID[v].str[9] = '0';
      m_MapID[v].length = 2;
      continue;
    }
    int idx = 11;
    m_MapID[v].str[--idx] = ',';
    while (x) {
      m_MapID[v].str[--idx] = x % 10 + '0';
      x /= 10;
    }
    m_MapID[v].start = idx;
    m_MapID[v].length = 11 - idx;
  }
}
void XJBG::SaveAnswer(int pid) {
#ifdef TEST
  std::cerr << m_debug[pid] << pid << ": SaveStart\n";
#endif
  int answers = *m_TotalAnswersPtr;
  FILE *fp;

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

    fp = fopen(RESULT, "w");
    fwrite(firBuf + firIdx, 1, 12 - firIdx, fp);
    fclose(fp);
    *m_FlagSavePtr = 0;
  }

  while (*m_FlagSavePtr != pid * 5) usleep(1);
  fp = fopen(RESULT, "at");
  fwrite(m_Answer0, 1, *m_AnswerLength0, fp);
  fclose(fp);
  if (pid == NTHREAD - 1) {
    *m_FlagSavePtr = 1;
  } else {
    *m_FlagSavePtr = (pid + 1) * 5;
  }

  while (*m_FlagSavePtr != pid * 5 + 1) usleep(1);
  fp = fopen(RESULT, "at");
  fwrite(m_Answer1, 1, *m_AnswerLength1, fp);
  fclose(fp);
  if (pid == NTHREAD - 1) {
    *m_FlagSavePtr = 2;
  } else {
    *m_FlagSavePtr = (pid + 1) * 5 + 1;
  }

  while (*m_FlagSavePtr != pid * 5 + 2) usleep(1);
  fp = fopen(RESULT, "at");
  fwrite(m_Answer2, 1, *m_AnswerLength2, fp);
  fclose(fp);
  if (pid == NTHREAD - 1) {
    *m_FlagSavePtr = 3;
  } else {
    *m_FlagSavePtr = (pid + 1) * 5 + 2;
  }

  while (*m_FlagSavePtr != pid * 5 + 3) usleep(1);
  fp = fopen(RESULT, "at");
  fwrite(m_Answer3, 1, *m_AnswerLength3, fp);
  fclose(fp);
  if (pid == NTHREAD - 1) {
    *m_FlagSavePtr = 4;
  } else {
    *m_FlagSavePtr = (pid + 1) * 5 + 3;
  }

  while (*m_FlagSavePtr != pid * 5 + 4) usleep(1);
  fp = fopen(RESULT, "at");
  fwrite(m_Answer4, 1, *m_AnswerLength4, fp);
  fclose(fp);
  if (pid == NTHREAD - 1) {
    *m_FlagSavePtr = 5;
  } else {
    *m_FlagSavePtr = (pid + 1) * 5 + 4;
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
    if (x == NTHREAD) {
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

#ifdef XJBGJUDGE
  delete xjbg;
#endif

  if (pid != 0) exit(0);

  return 0;
}