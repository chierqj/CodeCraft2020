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

struct PreBuffer {
  char str[10];
  uint32_t length;
};
const uint32_t MAXN = 50000 + 7;  // 总点数
const uint32_t NTHREAD = 8;       // 线程个数
const uint32_t ST[8] = {0, 6250, 9980, 10001, 11002, 17000, 25000, 28000};
const uint32_t ED[8] = {6250, 9980, 10001, 11002, 17000, 25000, 28000, MAXN};
uint32_t m_TotalAnswers = shmget(IPC_PRIVATE, 4, IPC_CREAT | 0600);
uint32_t m_FlagFork = shmget(IPC_PRIVATE, NTHREAD * 4, IPC_CREAT | 0600);
uint32_t m_FlagSave = shmget(IPC_PRIVATE, 4, IPC_CREAT | 0600);
uint32_t *m_TotalAnswersPtr;                       // share answers
uint32_t *m_FlagForkPtr;                           // share flag
uint32_t *m_FlagSavePtr;                           // share save
uint32_t m_maxID = 0;                              // 最大点
uint32_t m_answers = 0;                            // 环个数
uint32_t m_TempCount = 0;                          // 反向点的个数
uint32_t m_AnswerLength0 = 0;                      // 长度为3的环
uint32_t m_AnswerLength1 = 0;                      // 长度为4的环
uint32_t m_AnswerLength2 = 0;                      // 长度为5的环
uint32_t m_AnswerLength3 = 0;                      // 长度为6的环
uint32_t m_AnswerLength4 = 0;                      // 长度为7的环
uint32_t m_ThreeLength = 0;                        // 3环长度
char *m_ThreeCircle = new char[18];                // 长度为3的环
uint32_t m_TempPoint[2500];                        // 反向存那些点
char *m_Reachable = new char[MAXN];                // 可达
PreBuffer m_MapID[MAXN];                           // 预处理答案
uint32_t m_CountSons[MAXN], m_CountFathers[MAXN];  // 边数目
uint32_t m_Sons[MAXN][25];                         // 子结点
uint32_t m_Fathers[MAXN][25];                      // 父节点
char *m_Answer0 = new char[3000000 * 18];          // 长度为3的环
char *m_Answer1 = new char[3000000 * 24];          // 长度为4的环
char *m_Answer2 = new char[3000000 * 30];          // 长度为5的环
char *m_Answer3 = new char[3000000 * 36];          // 长度为6的环
char *m_Answer4 = new char[3000000 * 42];          // 长度为7的环

inline void GetShmPtr() {
  m_TotalAnswersPtr = (uint32_t *)shmat(m_TotalAnswers, NULL, 0);
  m_FlagForkPtr = (uint32_t *)shmat(m_FlagFork, NULL, 0);
  m_FlagSavePtr = (uint32_t *)shmat(m_FlagSave, NULL, 0);
}

void ParseInteger(uint32_t val) {
  uint32_t x = val;
  auto &mpid = m_MapID[val];
  if (x == 0) {
    mpid.str[0] = '0';
    mpid.str[1] = ',';
    mpid.length = 2;
    return;
  }
  char tmp[7];
  uint32_t idx = 7;
  tmp[--idx] = ',';
  while (x) {
    tmp[--idx] = x % 10 + '0';
    x /= 10;
  }
  memcpy(mpid.str, tmp + idx, 7 - idx);
  mpid.length = 7 - idx;
}

void handleLoadData(int start, int end) {
  for (uint32_t i = start; i < end; ++i) {
    if (m_CountSons[i] > 0 && m_CountFathers[i] > 0) {
      std::sort(m_Sons[i], m_Sons[i] + m_CountSons[i]);
      std::sort(m_Fathers[i], m_Fathers[i] + m_CountFathers[i],
                [&](const uint32_t &a, const uint32_t &b) { return a > b; });
      ParseInteger(i);
    }
  }
}

inline void addEdge(uint32_t u, uint32_t v) {
  if (u >= MAXN - 7 || v >= MAXN - 7) return;
  m_Sons[u][m_CountSons[u]++] = v;
  m_Fathers[v][m_CountFathers[v]++] = u;
  m_maxID = std::max(m_maxID, std::max(u, v) + 1);
}

void LoadData() {
  struct stat sb;
  uint32_t fd = open(TRAIN, O_RDONLY);
  fstat(fd, &sb);
  long long bufsize = sb.st_size;
  char *buffer = (char *)mmap(NULL, bufsize, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);

  uint32_t u = 0, v = 0, w = 0;
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
    while (*ptr != '\n') ++ptr;
    ++ptr;
    addEdge(u, v);
    u = v = 0;
  }

  std::thread Th[4];
  uint32_t st = 0, block = m_maxID / 4;
  for (int i = 0; i < 4; ++i) {
    int ed = (i == 3 ? m_maxID : st + block);
    Th[i] = std::thread(handleLoadData, st, ed);
    st += block;
  }
  for (int i = 0; i < 4; ++i) Th[i].join();

#ifdef LOCAL
  std::cerr << "@ LoadData: " << m_maxID << "\n";
#endif
}

void BackSearch(uint32_t st) {
  uint32_t i = 0, j = 0, k = 0;
  uint32_t count0 = m_CountFathers[st], count1, count2;
  uint32_t v1, v2, v3;
  uint32_t *father0 = m_Fathers[st], *father1, *father2;
  m_Reachable[st] = 7;
  m_TempPoint[m_TempCount++] = st;
  while (i++ < count0) {
    v1 = *(father0++);
    if (v1 <= st) break;
    m_Reachable[v1] |= 7;
    m_TempPoint[m_TempCount++] = v1;
    count1 = m_CountFathers[v1];
    father1 = m_Fathers[v1];
    j = 0;
    while (j++ < count1) {
      v2 = *(father1++);
      if (v2 <= st) break;
      m_Reachable[v2] |= 6;
      m_TempPoint[m_TempCount++] = v2;
      count2 = m_CountFathers[v2];
      father2 = m_Fathers[v2];
      k = 0;
      while (k++ < count2) {
        v3 = *(father2++);
        if (v3 <= st) break;
        m_Reachable[v3] |= 4;
        m_TempPoint[m_TempCount++] = v3;
      }
    }
  }
}

void doFindCircle(uint32_t st) {
  const auto &mpid0 = m_MapID[st];
  memcpy(m_ThreeCircle, mpid0.str, mpid0.length);
  uint32_t idx0 = mpid0.length, idx1 = 0;
  uint32_t it1 = 0, it2 = 0, it3 = 0, it4 = 0, it5 = 0, it6 = 0;
  uint32_t count0 = m_CountSons[st], count1, count2, count3, count4, count5;
  uint32_t v1, v2, v3, v4, v5, v6, v7;
  uint32_t *son0 = m_Sons[st], *son1, *son2, *son3, *son4, *son5;
  while (it1++ < count0) {
    v1 = *(son0++);
    if (v1 <= st) continue;
    const auto &mpid1 = m_MapID[v1];
    memcpy(m_ThreeCircle + idx0, mpid1.str, mpid1.length);
    idx1 = idx0 + mpid1.length;
    count1 = m_CountSons[v1];
    son1 = m_Sons[v1];
    it2 = 0;
    while (it2++ < count1) {
      v2 = *(son1++);
      if (v2 <= st) continue;
      const auto &mpid2 = m_MapID[v2];
      memcpy(m_ThreeCircle + idx1, mpid2.str, mpid2.length);
      m_ThreeLength = idx1 + mpid2.length;
      count2 = m_CountSons[v2];
      son2 = m_Sons[v2];
      it3 = 0;
      while (it3++ < count2) {
        v3 = *(son2++);
        if (v3 < st || v3 == v1) {
          continue;
        } else if (v3 == st) {
          memcpy(m_Answer0 + m_AnswerLength0, m_ThreeCircle, m_ThreeLength);
          m_AnswerLength0 += m_ThreeLength;
          *(m_Answer0 + m_AnswerLength0 - 1) = '\n';
          ++m_answers;
          continue;
        }
        count3 = m_CountSons[v3];
        son3 = m_Sons[v3];
        it4 = 0;
        const auto &mpid3 = m_MapID[v3];
        while (it4++ < count3) {
          v4 = *(son3++);
          if (!(m_Reachable[v4] & 4)) {
            continue;
          } else if (v4 == st) {
            memcpy(m_Answer1 + m_AnswerLength1, m_ThreeCircle, m_ThreeLength);
            m_AnswerLength1 += m_ThreeLength;
            memcpy(m_Answer1 + m_AnswerLength1, mpid3.str, mpid3.length);
            m_AnswerLength1 += mpid3.length;
            *(m_Answer1 + m_AnswerLength1 - 1) = '\n';
            ++m_answers;
            continue;
          } else if (v1 == v4 || v2 == v4) {
            continue;
          }
          count4 = m_CountSons[v4];
          son4 = m_Sons[v4];
          it5 = 0;
          const auto &mpid4 = m_MapID[v4];
          while (it5++ < count4) {
            v5 = *(son4++);
            if (!(m_Reachable[v5] & 2)) {
              continue;
            } else if (v5 == st) {
              memcpy(m_Answer2 + m_AnswerLength2, m_ThreeCircle, m_ThreeLength);
              m_AnswerLength2 += m_ThreeLength;
              memcpy(m_Answer2 + m_AnswerLength2, mpid3.str, mpid3.length);
              m_AnswerLength2 += mpid3.length;
              memcpy(m_Answer2 + m_AnswerLength2, mpid4.str, mpid4.length);
              m_AnswerLength2 += mpid4.length;
              *(m_Answer2 + m_AnswerLength2 - 1) = '\n';
              ++m_answers;
              continue;
            } else if (v1 == v5 || v2 == v5 || v3 == v5) {
              continue;
            }
            count5 = m_CountSons[v5];
            son5 = m_Sons[v5];
            it6 = 0;
            const auto &mpid5 = m_MapID[v5];
            while (it6++ < count5) {
              v6 = *(son5++);
              if (!(m_Reachable[v6] & 1)) {
                continue;
              } else if (v6 == st) {
                memcpy(m_Answer3 + m_AnswerLength3, m_ThreeCircle,
                       m_ThreeLength);
                m_AnswerLength3 += m_ThreeLength;
                memcpy(m_Answer3 + m_AnswerLength3, mpid3.str, mpid3.length);
                m_AnswerLength3 += mpid3.length;
                memcpy(m_Answer3 + m_AnswerLength3, mpid4.str, mpid4.length);
                m_AnswerLength3 += mpid4.length;
                memcpy(m_Answer3 + m_AnswerLength3, mpid5.str, mpid5.length);
                m_AnswerLength3 += mpid5.length;
                *(m_Answer3 + m_AnswerLength3 - 1) = '\n';
                ++m_answers;
                continue;
              } else if (v1 == v6 || v2 == v6 || v3 == v6 || v4 == v6) {
                continue;
              }
              const auto &mpid6 = m_MapID[v6];
              memcpy(m_Answer4 + m_AnswerLength4, m_ThreeCircle, m_ThreeLength);
              m_AnswerLength4 += m_ThreeLength;
              memcpy(m_Answer4 + m_AnswerLength4, mpid3.str, mpid3.length);
              m_AnswerLength4 += mpid3.length;
              memcpy(m_Answer4 + m_AnswerLength4, mpid4.str, mpid4.length);
              m_AnswerLength4 += mpid4.length;
              memcpy(m_Answer4 + m_AnswerLength4, mpid5.str, mpid5.length);
              m_AnswerLength4 += mpid5.length;
              memcpy(m_Answer4 + m_AnswerLength4, mpid6.str, mpid6.length);
              m_AnswerLength4 += mpid6.length;
              *(m_Answer4 + m_AnswerLength4 - 1) = '\n';
              ++m_answers;
            }
          }
        }
      }
    }
  }
}

void FindCircle(uint32_t pid) {
  uint32_t st = ST[pid];
  uint32_t ed = std::min(ED[pid], m_maxID);
  for (uint32_t i = st; i < ed; ++i) {
    if (m_CountSons[i] > 0 && m_CountFathers[i] > 0) {
      m_TempCount = 0;
      BackSearch(i);
      doFindCircle(i);
      for (uint32_t j = 0; j < m_TempCount; ++j) {
        uint32_t &p = m_TempPoint[j];
        m_Reachable[p] = 0;
      }
    }
  }
  *m_TotalAnswersPtr += m_answers;
  m_FlagForkPtr[pid] = 1;
}

void SaveAnswer(uint32_t pid) {
  uint32_t answers = *m_TotalAnswersPtr;
  FILE *fp;

  if (pid == 0) {
    uint32_t firIdx = 12;
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
  fwrite(m_Answer0, 1, m_AnswerLength0, fp);
  fclose(fp);
  if (pid == NTHREAD - 1) {
    *m_FlagSavePtr = 1;
  } else {
    *m_FlagSavePtr = (pid + 1) * 5;
  }

  while (*m_FlagSavePtr != pid * 5 + 1) usleep(1);
  fp = fopen(RESULT, "at");
  fwrite(m_Answer1, 1, m_AnswerLength1, fp);
  fclose(fp);
  if (pid == NTHREAD - 1) {
    *m_FlagSavePtr = 2;
  } else {
    *m_FlagSavePtr = (pid + 1) * 5 + 1;
  }

  while (*m_FlagSavePtr != pid * 5 + 2) usleep(1);
  fp = fopen(RESULT, "at");
  fwrite(m_Answer2, 1, m_AnswerLength2, fp);
  fclose(fp);
  if (pid == NTHREAD - 1) {
    *m_FlagSavePtr = 3;
  } else {
    *m_FlagSavePtr = (pid + 1) * 5 + 2;
  }

  while (*m_FlagSavePtr != pid * 5 + 3) usleep(1);
  fp = fopen(RESULT, "at");
  fwrite(m_Answer3, 1, m_AnswerLength3, fp);
  fclose(fp);
  if (pid == NTHREAD - 1) {
    *m_FlagSavePtr = 4;
  } else {
    *m_FlagSavePtr = (pid + 1) * 5 + 3;
  }

  while (*m_FlagSavePtr != pid * 5 + 4) usleep(1);
  fp = fopen(RESULT, "at");
  fwrite(m_Answer4, 1, m_AnswerLength4, fp);
  fclose(fp);
  if (pid == NTHREAD - 1) {
    *m_FlagSavePtr = 5;
  } else {
    *m_FlagSavePtr = (pid + 1) * 5 + 4;
  }

#ifdef LOCAL
  if (pid == 0) {
    std::cerr << "TotalAnswer: " << *m_TotalAnswersPtr << "\n";
  }
#endif
}

void WaitFork() {
  while (true) {
    uint32_t x = 0;
    for (uint32_t i = 0; i < NTHREAD; ++i) x += m_FlagForkPtr[i];
    if (x == NTHREAD) {
      return;
    }
    usleep(1);
  }
}

int main() {
  std::cerr << std::fixed << std::setprecision(3);

  LoadData();

  pid_t Children[NTHREAD - 1] = {0};
  uint32_t pid = 0;
  for (uint32_t i = 0; i < NTHREAD - 1; i++) {
    if (pid == 0) {
      uint32_t x = fork();
      Children[i] = x;
      if (x == 0) {
        pid = i + 1;
        break;
      }
    }
  }

  GetShmPtr();
  FindCircle(pid);
  WaitFork();
  SaveAnswer(pid);

#ifdef XJBGJUDGE
  delete[] m_Answer0;
  delete[] m_Answer1;
  delete[] m_Answer2;
  delete[] m_Answer3;
  delete[] m_Answer4;
  delete[] m_ThreeCircle;
  delete[] m_Reachable;
#endif

  if (pid != 0) exit(0);
  exit(0);
}