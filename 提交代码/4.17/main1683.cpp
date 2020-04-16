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

class XJBG {
 public:
#ifdef XJBGJUDGE
  ~XJBG() {
    delete[] m_Answer0;
    delete[] m_Answer1;
    delete[] m_Answer2;
    delete[] m_Answer3;
    delete[] m_Answer4;
    delete[] m_ThreeCircle;
    delete[] m_Reachable;
    delete m_AnswerLength0;
    delete m_AnswerLength1;
    delete m_AnswerLength2;
    delete m_AnswerLength3;
    delete m_AnswerLength4;
  }
#endif

 private:
  struct PreBuffer {
    char str[10];
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
  void WaitFork();           // 等待进程
  void GetShmPtr();          // 获取共享数据句柄

 private:
  inline void addEdge(int u, int v, int w);  // 加边
  void BackSearch(int st);                   // 反向剪枝
  void doFindCircle(int st);                 // 正向找环
  inline void spliceMid(char *ans, uint32_t *len, int node);
  inline void spliceEnd(char *ans, uint32_t *len, int node);
  inline void spliceBegin(char *ans, uint32_t *cur, uint32_t len);
  inline void spliceBeginEnd(char *ans, uint32_t *cur, uint32_t len);

 public:
  static const int MAXN = 200000 + 7;  // 总点数
  static const int NTHREAD = 16;       // 线程个数
  const int PARAM[NTHREAD] = {1, 2,  3,  4,  5,  6,  7,  8,
                              9, 10, 11, 12, 13, 14, 15, 16};  // 线程权重
  const int PARAM_SUM = 136;

 private:
  int m_maxID = 0;                              // 最大点
  int m_Sons[MAXN][30];                         // 子结点
  int m_Fathers[MAXN][30];                      // 父节点
  int m_CountSons[MAXN], m_CountFathers[MAXN];  // 边数目
  int m_edgeNum = 0;                            // 边数目
  PreBuffer m_MapID[MAXN];                      // 预处理答案
  int m_Jobs[MAXN];                             // 点集合
  int m_JobsCount = 0;
  int m_answers = 0;                         // 环个数
  char *m_Reachable;                         // 可达
  char *m_Answer0 = new char[3000000 * 21];  // 长度为3的环
  char *m_Answer1 = new char[3000000 * 28];  // 长度为4的环
  char *m_Answer2 = new char[3000000 * 35];  // 长度为5的环
  char *m_Answer3 = new char[3000000 * 42];  // 长度为6的环
  char *m_Answer4 = new char[3000000 * 49];  // 长度为7的环
  uint32_t *m_AnswerLength0 = new uint32_t;  // 长度为3的环buffer长度
  uint32_t *m_AnswerLength1 = new uint32_t;  // 长度为4的环buffer长度
  uint32_t *m_AnswerLength2 = new uint32_t;  // 长度为5的环buffer长度
  uint32_t *m_AnswerLength3 = new uint32_t;  // 长度为6的环buffer长度
  uint32_t *m_AnswerLength4 = new uint32_t;  // 长度为7的环buffer长度
  char *m_ThreeCircle = new char[22];        // 前缀长度为3的环
  int m_TempPoint[3000];                     // 反向存那些点
  int m_TempCount = 0;                       // 反向点的个数

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
      // w = w * 10 + *ptr - '0';
      ++ptr;
    }
    ++ptr;
    addEdge(u, v, w);
    u = v = w = 0;
  }

  for (int i = 0; i < m_maxID; ++i) {
    if (m_CountSons[i] > 0 && m_CountFathers[i] > 0) {
      std::sort(m_Sons[i], m_Sons[i] + m_CountSons[i]);
      m_Jobs[m_JobsCount++] = i;
    }
  }
  m_Reachable = new char[m_maxID];

#ifdef TEST
  std::cerr << "@ LoadData: (V: " << m_maxID << ", E: " << m_edgeNum << ") #";
  t.LogTime();
#endif
}

void XJBG::BackSearch(int st) {
  int v1, v2, v3;
  int count0 = m_CountFathers[st], count1, count2;
  int *father0 = m_Fathers[st], *father1, *father2;
  m_Reachable[st] = 7;
  m_TempPoint[m_TempCount++] = st;

  for (int i = 0; i < count0; ++i) {
    v1 = *(father0 + i);
    if (v1 <= st) continue;
    m_Reachable[v1] |= 7;
    m_TempPoint[m_TempCount++] = v1;
    count1 = m_CountFathers[v1];
    father1 = m_Fathers[v1];
    for (int j = 0; j < count1; ++j) {
      v2 = *(father1 + j);
      if (v2 <= st) continue;
      m_Reachable[v2] |= 6;
      m_TempPoint[m_TempCount++] = v2;
      count2 = m_CountFathers[v2];
      father2 = m_Fathers[v2];
      for (int k = 0; k < count2; ++k) {
        v3 = *(father2 + k);
        if (v3 <= st || v3 == v1) continue;
        m_Reachable[v3] |= 4;
        m_TempPoint[m_TempCount++] = v3;
      }
    }
  }
}

// 拼接开始保存的三个点
inline void XJBG::spliceBegin(char *ans, uint32_t *cur, uint32_t len) {
  memcpy(ans + *cur, m_ThreeCircle, len);
  *cur += len;
}
// 恰好保存的三个点是一个三环
inline void XJBG::spliceBeginEnd(char *ans, uint32_t *cur, uint32_t len) {
  memcpy(ans + *cur, m_ThreeCircle, len);
  *cur += len;
  *(ans + *cur - 1) = '\n';
  ++m_answers;
}
// 拼接中间结点
inline void XJBG::spliceMid(char *ans, uint32_t *len, int node) {
  const auto &mpid = m_MapID[node];
  memcpy(ans + *len, mpid.str, mpid.length);
  *len += mpid.length;
}
// 拼接最后一个点，注意换行符
inline void XJBG::spliceEnd(char *ans, uint32_t *len, int node) {
  const auto &mpid = m_MapID[node];
  memcpy(ans + *len, mpid.str, mpid.length);
  *len += mpid.length;
  *(ans + *len - 1) = '\n';
  ++m_answers;
}

void XJBG::doFindCircle(int st) {
  const auto &mpid0 = m_MapID[st];
  memcpy(m_ThreeCircle, mpid0.str, mpid0.length);
  int IDX0 = mpid0.length, IDX1 = 0, IDX2 = 0;

  int v1, v2, v3, v4, v5, v6, v7;
  int count0 = m_CountSons[st], count1, count2, count3, count4, count5;
  int *son0 = m_Sons[st], *son1, *son2, *son3, *son4, *son5;

  for (int it1 = 0; it1 < count0; ++it1) {
    v1 = *(son0 + it1);
    if (v1 < st) continue;

    const auto &mpid1 = m_MapID[v1];
    memcpy(m_ThreeCircle + IDX0, mpid1.str, mpid1.length);
    IDX1 = IDX0 + mpid1.length;

    count1 = m_CountSons[v1];
    son1 = m_Sons[v1];
    for (int it2 = 0; it2 < count1; ++it2) {
      v2 = *(son1 + it2);
      if (v2 <= st) continue;

      const auto &mpid2 = m_MapID[v2];
      memcpy(m_ThreeCircle + IDX1, mpid2.str, mpid2.length);
      IDX2 = IDX1 + mpid2.length;

      count2 = m_CountSons[v2];
      son2 = m_Sons[v2];
      for (int it3 = 0; it3 < count2; ++it3) {
        v3 = *(son2 + it3);
        if (v3 < st) {
          continue;
        } else if (!(v3 ^ st)) {
          spliceBeginEnd(m_Answer0, m_AnswerLength0, IDX2);
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
            spliceBegin(m_Answer1, m_AnswerLength1, IDX2);
            spliceEnd(m_Answer1, m_AnswerLength1, v3);
            continue;
          } else if (!(v4 ^ v1) || !(v4 ^ v2)) {
            continue;
          }
          count4 = m_CountSons[v4];
          son4 = m_Sons[v4];
          for (int it5 = 0; it5 < count4; ++it5) {
            v5 = *(son4 + it5);
            if (!(m_Reachable[v5] & 2)) {
              continue;
            } else if (!(v5 ^ st)) {
              spliceBegin(m_Answer2, m_AnswerLength2, IDX2);
              spliceMid(m_Answer2, m_AnswerLength2, v3);
              spliceEnd(m_Answer2, m_AnswerLength2, v4);
              continue;
            } else if (!(v5 ^ v1) || !(v5 ^ v2) || !(v5 ^ v3)) {
              continue;
            }
            count5 = m_CountSons[v5];
            son5 = m_Sons[v5];
            for (int it6 = 0; it6 < count5; ++it6) {
              v6 = *(son5 + it6);
              if (!(m_Reachable[v6] & 1)) {
                continue;
              } else if (!(v6 ^ st)) {
                spliceBegin(m_Answer3, m_AnswerLength3, IDX2);
                spliceMid(m_Answer3, m_AnswerLength3, v3);
                spliceMid(m_Answer3, m_AnswerLength3, v4);
                spliceEnd(m_Answer3, m_AnswerLength3, v5);
                continue;
              } else if (!(v6 ^ v1) || !(v6 ^ v2) || !(v6 ^ v3) || !(v6 ^ v4)) {
                continue;
              }
              spliceBegin(m_Answer4, m_AnswerLength4, IDX2);
              spliceMid(m_Answer4, m_AnswerLength4, v3);
              spliceMid(m_Answer4, m_AnswerLength4, v4);
              spliceMid(m_Answer4, m_AnswerLength4, v5);
              spliceEnd(m_Answer4, m_AnswerLength4, v6);
            }
          }
        }
      }
    }
  }
}

void XJBG::FindCircle(int pid) {
  ScopeTime t;

#ifdef TEST
  std::cerr << m_debug[pid] << pid << ": FindStart\n";
#endif
  int block = m_JobsCount / PARAM_SUM, st = 0;
  for (int i = 0; i < pid; ++i) st += PARAM[i] * block;
  int ed = (pid == NTHREAD - 1 ? m_JobsCount : st + block * PARAM[pid]);

  for (int i = st; i < ed; ++i) {
    int v = m_Jobs[i];
    m_TempCount = 0;
    this->BackSearch(v);
    this->doFindCircle(v);
    for (int j = 0; j < m_TempCount; ++j) {
      m_Reachable[m_TempPoint[j]] = 0;
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
  // m_MapID.resize(m_maxID);

  for (int i = 0; i < m_JobsCount; ++i) {
    int v = m_Jobs[i];
    int x = v;
    auto &mpid = m_MapID[v];
    if (x == 0) {
      char tmp[2] = {'0', ','};
      memcpy(mpid.str, tmp, 2);
      mpid.length = 2;
      continue;
    }
    char tmp[10];
    int idx = 10;
    tmp[--idx] = ',';
    while (x) {
      tmp[--idx] = x % 10 + '0';
      x /= 10;
    }
    memcpy(mpid.str, tmp + idx, 10 - idx);
    mpid.length = 10 - idx;
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