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
#include <mutex>
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
    for (int i = 0; i < NTHREAD; ++i) {
      delete[] ThDatas[i].Reachable;
    }
  }
#endif

 public:
  static const int MAXN = 200000 + 7;  // 总点数
  static const int NTHREAD = 4;        // 线程个数

 private:
  // int转char*
  struct PreBuffer {
    char *str;
    char *str1;
    int length;
  };
  // 线程数据
  struct ThreadData {
    int NowJob = 0;                    // 当前任务id
    int answers = 0;                   // 环数
    char *Reachable = new char[MAXN];  // 反向标记
    int TempPoint[3000];               // 反向修改的点
    int TempPointCount = 0;            // 反向修改的点的个数
    int DFSCOUNT = 0;                  // 递归次数
  };
  // 按照起点保存答案
  struct Circle {
    char **buffer[5];
    int cur[5];
    int count[5];
  };

 public:
  void LoadAndPreTreatData();  // 加载数据
  void FindCircle();           // 找环
  void SaveAnswer();           // 保存答案

 private:
  inline void addEdge(int u, int v, int w);  // 加边
  void parseIntger(int v);                   // 解析整数
  void handleThread(int pid);                // 处理线程
  void getNextJobID(int &job);               // 获取任务id

  void BackSearch(ThreadData &Data);    // 反向剪枝
  void doFindCircle(ThreadData &Data);  // 正向找环
  inline void rellocBuffer(Circle &circle, int dep);
  inline void spliceNode(Circle &circle, int dep, const int &node);
  inline void spliceEnd(Circle &circle, int dep, const int &node);

 private:
  int m_maxID = 0;                              // 最大点
  int m_Sons[MAXN][30];                         // 子结点
  int m_Fathers[MAXN][30];                      // 父节点
  int m_CountSons[MAXN], m_CountFathers[MAXN];  // 边数目
  int m_edgeNum = 0;                            // 边数目
  PreBuffer m_MapID[MAXN];                      // 预处理答案
  int m_Jobs[MAXN];                             // 点集合
  int m_JobCount = 0;                           // 点数目
  int m_JobCur = 0;                             // 当前处理到第几个job
  std::mutex m_mtx;                             // 线程锁
  ThreadData ThDatas[NTHREAD];                  // 线程数据
  int m_answers = 0;                            // 总答案
  Circle AllCircles[MAXN];                      // 按照起点保存的环
};

inline void XJBG::addEdge(int u, int v, int w) {
  m_Sons[u][m_CountSons[u]++] = v;
  m_Fathers[v][m_CountFathers[v]++] = u;
  ++m_edgeNum;
  m_maxID = std::max(m_maxID, std::max(u, v) + 1);
}
void XJBG::parseIntger(int x) {
  auto &mpid = m_MapID[x];
  if (x == 0) {
    mpid.str = new char[2];
    mpid.str1 = new char[2];
    mpid.str1[0] = '0';
    mpid.str1[1] = ',';
    mpid.str1[0] = '0';
    mpid.str1[1] = '\n';
    mpid.length = 2;
    return;
  }
  int val = x, idx = 10;
  char tmp[10];
  tmp[--idx] = ',';
  while (val) {
    tmp[--idx] = val % 10 + '0';
    val /= 10;
  }
  mpid.str = new char[10 - idx];
  mpid.str1 = new char[10 - idx];
  memcpy(mpid.str, tmp + idx, 10 - idx);
  tmp[9] = '\n';
  memcpy(mpid.str1, tmp + idx, 10 - idx);
  mpid.length = 10 - idx;
}
void XJBG::LoadAndPreTreatData() {
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
      m_Jobs[m_JobCount++] = i;
      this->parseIntger(i);
    }
  }

#ifdef TEST
  std::cerr << "@ LoadData: (V: " << m_maxID << ", E: " << m_edgeNum
            << ", Job: " << m_JobCount << ") #";
  t.LogTime();
#endif
}

void XJBG::BackSearch(ThreadData &Data) {
  int v1, v2, v3;
  int count0 = m_CountFathers[Data.NowJob], count1, count2;
  int *father0 = m_Fathers[Data.NowJob], *father1, *father2;
  Data.Reachable[Data.NowJob] = 7;
  Data.TempPoint[Data.TempPointCount++] = Data.NowJob;
  for (int i = 0; i < count0; ++i) {
    v1 = *(father0 + i);
    if (v1 <= Data.NowJob) continue;
    Data.Reachable[v1] |= 7;
    Data.TempPoint[Data.TempPointCount++] = v1;
    count1 = m_CountFathers[v1];
    father1 = m_Fathers[v1];
    for (int j = 0; j < count1; ++j) {
      v2 = *(father1 + j);
      if (v2 <= Data.NowJob) continue;
      Data.Reachable[v2] |= 6;
      Data.TempPoint[Data.TempPointCount++] = v2;
      count2 = m_CountFathers[v2];
      father2 = m_Fathers[v2];
      for (int k = 0; k < count2; ++k) {
        v3 = *(father2 + k);
        if (v3 <= Data.NowJob || v3 == v1) continue;
        Data.Reachable[v3] |= 4;
        Data.TempPoint[Data.TempPointCount++] = v3;
      }
    }
  }
}

inline void XJBG::rellocBuffer(Circle &circle, int dep) {
  if ((circle.count[dep] & -circle.count[dep]) == circle.count[dep]) {
    circle.buffer[dep] = (char **)realloc(
        circle.buffer[dep], (circle.count[dep] * 2 + 1) * sizeof(char **));
  }
}
inline void XJBG::spliceNode(Circle &circle, int dep, const int &node) {
  rellocBuffer(circle, dep);
  auto &mpid = m_MapID[node];
  circle.buffer[dep][circle.count[dep]++] = mpid.str;
}
inline void XJBG::spliceEnd(Circle &circle, int dep, const int &node) {
  rellocBuffer(circle, dep);
  auto &mpid = m_MapID[node];
  circle.buffer[dep][circle.count[dep]++] = mpid.str1;
}
void XJBG::doFindCircle(ThreadData &Data) {
  const auto &mpid0 = m_MapID[Data.NowJob];
  int v1, v2, v3, v4, v5, v6, v7;
  int count0 = m_CountSons[Data.NowJob], count1, count2, count3, count4, count5;
  int *son0 = m_Sons[Data.NowJob], *son1, *son2, *son3, *son4, *son5;

  auto &circle = AllCircles[Data.NowJob];

  for (int it1 = 0; it1 < count0; ++it1) {
    v1 = *(son0 + it1);
    if (v1 < Data.NowJob) continue;
    count1 = m_CountSons[v1];
    son1 = m_Sons[v1];
    for (int it2 = 0; it2 < count1; ++it2) {
      v2 = *(son1 + it2);
      if (v2 <= Data.NowJob) continue;
      count2 = m_CountSons[v2];
      son2 = m_Sons[v2];
      for (int it3 = 0; it3 < count2; ++it3) {
        v3 = *(son2 + it3);
        if (v3 < Data.NowJob) {
          continue;
        } else if (v3 == Data.NowJob) {
          spliceNode(circle, 0, Data.NowJob);
          spliceNode(circle, 0, v1);
          spliceEnd(circle, 0, v2);
          ++Data.answers;
          continue;
        } else if (v3 == v1) {
          continue;
        }
        count3 = m_CountSons[v3];
        son3 = m_Sons[v3];
        for (int it4 = 0; it4 < count3; ++it4) {
          v4 = *(son3 + it4);

          if (!(Data.Reachable[v4] & 4)) {
            continue;
          } else if (v4 == Data.NowJob) {
            spliceNode(circle, 1, Data.NowJob);
            spliceNode(circle, 1, v1);
            spliceNode(circle, 1, v2);
            spliceEnd(circle, 1, v3);
            ++Data.answers;
            continue;
          } else if (v4 == v1 || v4 == v2) {
            continue;
          }
          count4 = m_CountSons[v4];
          son4 = m_Sons[v4];
          for (int it5 = 0; it5 < count4; ++it5) {
            v5 = *(son4 + it5);

            if (!(Data.Reachable[v5] & 2)) {
              continue;
            } else if (v5 == Data.NowJob) {
              spliceNode(circle, 2, Data.NowJob);
              spliceNode(circle, 2, v1);
              spliceNode(circle, 2, v2);
              spliceNode(circle, 2, v3);
              spliceEnd(circle, 2, v4);
              ++Data.answers;
              continue;
            } else if (v5 == v1 || v5 == v2 || v5 == v3) {
              continue;
            }
            count5 = m_CountSons[v5];
            son5 = m_Sons[v5];
            for (int it6 = 0; it6 < count5; ++it6) {
              v6 = *(son5 + it6);
              if (!(Data.Reachable[v6] & 1)) {
                continue;
              } else if (v6 == Data.NowJob) {
                spliceNode(circle, 3, Data.NowJob);
                spliceNode(circle, 3, v1);
                spliceNode(circle, 3, v2);
                spliceNode(circle, 3, v3);
                spliceNode(circle, 3, v4);
                spliceEnd(circle, 3, v5);
                ++Data.answers;
                continue;
              } else if (v6 == v1 || v6 == v2 || v6 == v3 || v6 == v4 ||
                         v6 == v5) {
                continue;
              }
              spliceNode(circle, 4, Data.NowJob);
              spliceNode(circle, 4, v1);
              spliceNode(circle, 4, v2);
              spliceNode(circle, 4, v3);
              spliceNode(circle, 4, v4);
              spliceNode(circle, 4, v5);
              spliceEnd(circle, 4, v6);
              ++Data.answers;
            }
          }
        }
      }
    }
  }
}

void XJBG::getNextJobID(int &job) {
  m_mtx.lock();
  if (m_JobCur < m_JobCount) {
    job = m_Jobs[m_JobCur++];
  } else {
    job = -1;
  }
  m_mtx.unlock();
}
void XJBG::handleThread(int pid) {
  auto &Data = ThDatas[pid];
  while (true) {
    getNextJobID(Data.NowJob);
    if (Data.NowJob == -1) break;
    Data.TempPointCount = 0;
    this->BackSearch(Data);
    this->doFindCircle(Data);
    for (int j = 0; j < Data.TempPointCount; ++j) {
      Data.Reachable[Data.TempPoint[j]] = 0;
    }
  }
}
void XJBG::FindCircle() {
  std::thread Th[NTHREAD];
  for (int i = 0; i < NTHREAD; ++i) {
    Th[i] = std::thread(&XJBG::handleThread, this, i);
  }
  for (int i = 0; i < NTHREAD; ++i) Th[i].join();
  for (int i = 0; i < NTHREAD; ++i) m_answers += ThDatas[i].answers;
#ifdef LOCAL
  std::cerr << "@ circles: " << m_answers << "\n";
#endif
}

void XJBG::SaveAnswer() {
  return;
  int firIdx = 12, answers = m_answers;
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
  // fclose(fp);
}

int main() {
  std::cerr << std::fixed << std::setprecision(3);

  XJBG *xjbg = new XJBG();

  xjbg->LoadAndPreTreatData();
  xjbg->FindCircle();
  xjbg->SaveAnswer();

#ifdef XJBGJUDGE
  delete xjbg;
#endif
  return 0;
}