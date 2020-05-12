#include <fcntl.h>
#include <sys/mman.h>
#include <sys/time.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <numeric>
#include <queue>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#define uint uint32_t
#define W3_MAX 715827882
#define W5_MAX 429496729
#define P3(x) ((x << 1) + x)
#define P5(x) ((x << 2) + x)
#define P10(x) ((x << 3) + (x << 1))

#ifdef LOCAL
#define TRAIN "../data/18908526/test_data.txt"
#define RESULT "/dev/shm/result.txt"
#else
#define TRAIN "/data/test_data.txt"
#define RESULT "/projects/student/result.txt"
#endif
typedef std::pair<uint, uint> Pair;

/*
 * 常量定义
 */
const uint MAXEDGE = 2000000 + 7;  // 最多边数目
const uint MAXN = MAXEDGE << 1;    // 最多点数目
const uint NTHREAD = 4;            // 线程个数
const uint NUMLENGTH = 12;         // ID最大长度
const uint BUFFERBLOCK = 128;

struct DFSEdge {
  uint idx, w;
};
struct Edge {
  uint u, v, w;
};
struct PreBuffer {
  uint len;
  char str[NUMLENGTH];
};
struct Answer {
  std::vector<uint> cycle[5];
};
struct ThData {
  uint ReachCount = 0;    // 反向可达点数目
  char Reach[MAXN];       // 标记反向可达
  uint ReachPoint[MAXN];  // 可达点集合
  uint LastWeight[MAXN];  // 最后一步权重
};

uint MaxID = 0;                                     // 最大点
uint Answers = 0;                                   // 环个数
uint EdgesCount = 0;                                // 边数目
uint JobsCount = 0;                                 // 有效点数目
uint TotalBufferSize = 0;                           // 总buffer大小
uint FirstBufLen = 0;                               // 换个数bufsize
uint JobCur = 0;                                    // Job光标
uint OffSetCur = 0;                                 // OffSet光标
bool Done[BUFFERBLOCK];                             // 构造buffer
uint AnswerBufferLen[BUFFERBLOCK];                  // bufferlen
uint Head[MAXN], HeadLen[MAXN];                     // 前向图
uint Back[MAXN], BackLen[MAXN];                     // 后向图
DFSEdge G[MAXN];                                    // 前向图
DFSEdge GBack[MAXN];                                // 后向图
std::atomic_flag _JOB_LOCK_ = ATOMIC_FLAG_INIT;     // job lock
std::atomic_flag _OFFSET_LOCK_ = ATOMIC_FLAG_INIT;  // offset lock
char FirstBuf[NUMLENGTH];                           // 环个数buf
std::vector<std::array<uint, 4>> OffSet;  // 分块输出，每一块位置
uint Jobs[MAXN];                          // 有效点
PreBuffer MapID[MAXN];                    // 解析int
std::vector<Answer> Cycles;               // 所有答案
Edge Edges[MAXEDGE];                      // 读入正向边集合
ThData ThreadData[NTHREAD];               // 线程找环

char *AnswerBuffer[BUFFERBLOCK];

struct HashTable {
  static const int MOD1 = 6893911;
  static const int MOD2 = 5170427;

  struct Data {
    uint key;
    int val = -1;
  };

  uint count = 0;
  Data Map[MOD1];
  uint HashIdx[MOD1];

  uint size() { return count; }
  uint hash(const uint &k, const int &i) {
    return (k % MOD1 + i * (MOD2 - k % MOD2)) % MOD1;
  }
  void Insert(const uint &key) {
    for (int i = 0; i < MOD1; ++i) {
      const uint &val = this->hash(key, i);
      if (Map[val].val == -1) {
        Map[val].key = key;
        Map[val].val = MaxID++;
        HashIdx[count++] = val;
        break;
      }
      if (Map[val].key == key) break;
    }
  }
  int Query(const uint &key) {
    for (int i = 0; i < MOD1; ++i) {
      const uint &val = hash(key, i);
      if (Map[val].val == -1) return -1;
      if (Map[val].key == key) return Map[val].val;
    }
    return -1;
  }
  void Sort() {
    std::sort(HashIdx, HashIdx + count, [&](const uint &x, const uint &y) {
      return Map[x].key < Map[y].key;
    });
    for (int i = 0; i < count; ++i) {
      const uint &idx = HashIdx[i];
      Map[idx].val = i;
      auto &mpid = MapID[i];
      sprintf(mpid.str, "%d,", Map[idx].key);
      mpid.len = strlen(mpid.str);
    }
  }
};
HashTable HashMap;

inline void addEdge(const uint &u, const uint &v, const uint &w) {
  auto &e = Edges[EdgesCount++];
  e.u = u;
  e.v = v;
  e.w = w;
}

void CreateHashTable() {
#ifdef LOCAL
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif

  uint pre = 0;
  for (uint i = 0; i < EdgesCount; ++i) {
    const auto &e = Edges[i];
    HashMap.Insert(e.u);
  }
  HashMap.Sort();

#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ HashMap:\t[cost: %.4fs]\n", t4 - t1);
#endif
}

void CreateForwardGraph() {
#ifdef LOCAL
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif

  for (uint i = 0; i < EdgesCount; ++i) {
    const auto &e = Edges[i];
    uint v = HashMap.Query(e.v);
    if (v != -1) {
      uint u = HashMap.Query(e.u);
      ++HeadLen[u];
    }
  }
  Head[0] = 0;
  for (uint i = 1; i <= MaxID; ++i) Head[i] = Head[i - 1] + HeadLen[i - 1];
  std::vector<uint> cnt(MaxID, 0);
  for (uint i = 0; i < EdgesCount; ++i) {
    const auto &e = Edges[i];
    uint v = HashMap.Query(e.v);
    if (v != -1) {
      uint u = HashMap.Query(e.u);
      G[Head[u] + cnt[u]].idx = v;
      G[Head[u] + cnt[u]++].w = e.w;
    }
  }
#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ GHead:\t[cost: %.4fs]\n", t4 - t1);
#endif
}
void CreateBackGraph() {
#ifdef LOCAL
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif
  for (uint i = 0; i < EdgesCount; ++i) {
    const auto &e = Edges[i];
    uint v = HashMap.Query(e.v);
    if (v != -1) {
      ++BackLen[v];
    }
  }
  Back[0] = 0;
  for (uint i = 1; i <= MaxID; ++i) Back[i] = Back[i - 1] + BackLen[i - 1];
  std::vector<uint> cnt(MaxID, 0);
  for (uint i = 0; i < EdgesCount; ++i) {
    const auto &e = Edges[i];
    uint v = HashMap.Query(e.v);
    if (v != -1) {
      uint u = HashMap.Query(e.u);
      GBack[Back[v] + cnt[v]].idx = u;
      GBack[Back[v] + cnt[v]++].w = e.w;
    }
  }
#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ GBack:\t[cost: %.4fs]\n", t4 - t1);
#endif
}
void SortEdge(uint pid) {
  for (uint i = pid; i < MaxID; i += NTHREAD) {
    std::sort(
        G + Head[i], G + Head[i] + HeadLen[i],
        [&](const DFSEdge &e1, const DFSEdge &e2) { return e1.idx < e2.idx; });
    std::sort(
        GBack + Back[i], GBack + Back[i] + BackLen[i],
        [&](const DFSEdge &e1, const DFSEdge &e2) { return e1.idx > e2.idx; });
  }
}

void LoadData() {
#ifdef LOCAL
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif

  uint fd = open(TRAIN, O_RDONLY);
  uint bufsize = lseek(fd, 0, SEEK_END);
  char *buffer = (char *)mmap(NULL, bufsize, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);
  const char *ptr = buffer;
  const char *end = buffer + bufsize;

  uint u = 0, v = 0, w = 0;
  while (ptr < end) {
    while (*ptr != ',') {
      u = P10(u) + *ptr - '0';
      ++ptr;
    }
    ++ptr;
    while (*ptr != ',') {
      v = P10(v) + *ptr - '0';
      ++ptr;
    }
    ++ptr;
    while (*ptr != '\r' && *ptr != '\n') {
      w = P10(w) + *ptr - '0';
      ++ptr;
    }
    if (*ptr == '\r') ++ptr;
    ++ptr;
    addEdge(u, v, w);
    u = v = w = 0;
  }

#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ Buffer:\t[cost: %.4fs]\n", t4 - t1);
#endif

  CreateHashTable();
  std::thread Th[4];
  Th[0] = std::thread(CreateForwardGraph);
  Th[1] = std::thread(CreateBackGraph);
  Th[0].join();
  Th[1].join();
  for (uint i = 0; i < NTHREAD; ++i) Th[i] = std::thread(SortEdge, i);
  for (uint i = 0; i < NTHREAD; ++i) Th[i].join();

  for (uint i = 0; i < MaxID; ++i) {
    if (HeadLen[i] > 0 && BackLen[i] > 0) {
      Jobs[JobsCount++] = i;
    }
  }

#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t43 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ [U: %d, E: %d, Job: %d] [cost: %.4fs]\n\n", MaxID, EdgesCount,
         JobsCount, t43 - t1);
#endif
}

inline bool judge(const uint &w1, const uint &w2) {
  return (w2 > W5_MAX || w1 <= P5(w2)) && (w1 > W3_MAX || P3(w1) >= w2);
}

void BackSearch(ThData &Data, const uint &st) {
  for (uint i = 0; i < Data.ReachCount; ++i) {
    const uint &v = Data.ReachPoint[i];
    Data.Reach[v] = 0;
  }
  Data.ReachCount = 0;
  Data.ReachPoint[Data.ReachCount++] = st;
  Data.Reach[st] = 7;

  const DFSEdge *e1 = &GBack[Back[st]];
  for (uint it1 = Back[st]; it1 < Back[st + 1]; ++it1, ++e1) {
    const uint &v1 = e1->idx;
    if (v1 <= st) break;
    const uint &w1 = e1->w;
    Data.LastWeight[v1] = w1;
    Data.Reach[v1] = 7;
    Data.ReachPoint[Data.ReachCount++] = v1;

    const DFSEdge *e2 = &GBack[Back[v1]];
    for (uint it2 = Back[v1]; it2 < Back[v1 + 1]; ++it2, ++e2) {
      const uint &v2 = e2->idx;
      if (v2 <= st) break;
      const uint &w2 = e2->w;

      if (!judge(w2, w1)) continue;
      Data.Reach[v2] |= 6;
      Data.ReachPoint[Data.ReachCount++] = v2;

      const DFSEdge *e3 = &GBack[Back[v2]];
      for (uint it3 = Back[v2]; it3 < Back[v2 + 1]; ++it3, ++e3) {
        const uint &v3 = e3->idx;
        if (v3 <= st) break;
        if (v3 == v1) continue;
        const uint &w3 = e3->w;
        if (!judge(w3, w2)) continue;
        Data.Reach[v3] |= 4;
        Data.ReachPoint[Data.ReachCount++] = v3;
      }
    }
  }
}

void ForwardSearch(ThData &Data, const uint &st) {
  std::vector<uint>(&ret)[5] = Cycles[st].cycle;

  const DFSEdge *e1 = &G[Head[st]];
  for (uint it1 = Head[st]; it1 < Head[st + 1]; ++it1, ++e1) {
    const uint &v1 = e1->idx, &w1 = e1->w;
    if (v1 < st) continue;

    const DFSEdge *e2 = &G[Head[v1]];

    for (uint it2 = Head[v1]; it2 < Head[v1 + 1]; ++it2, ++e2) {
      const uint &v2 = e2->idx, &w2 = e2->w;
      if (v2 <= st || !judge(w1, w2)) continue;

      const DFSEdge *e3 = &G[Head[v2]];

      for (uint it3 = Head[v2]; it3 < Head[v2 + 1]; ++it3, ++e3) {
        const uint &v3 = e3->idx, &w3 = e3->w;

        if (v3 < st || v3 == v1 || !judge(w2, w3)) {
          continue;
        } else if (v3 == st) {
          if (judge(w3, w1)) {
            ret[0].insert(ret[0].end(), {st, v1, v2});
          }
          continue;
        }

        const DFSEdge *e4 = &G[Head[v3]];

        for (uint it4 = Head[v3]; it4 < Head[v3 + 1]; ++it4, ++e4) {
          const uint &v4 = e4->idx, &w4 = e4->w;

          if (!(Data.Reach[v4] & 4) || v1 == v4 || v2 == v4) {
            continue;
          } else if (!judge(w3, w4)) {
            continue;
          } else if (v4 == st) {
            if (judge(w4, w1)) {
              ret[1].insert(ret[1].end(), {st, v1, v2, v3});
            }
            continue;
          }

          const DFSEdge *e5 = &G[Head[v4]];

          for (uint it5 = Head[v4]; it5 < Head[v4 + 1]; ++it5, ++e5) {
            const uint &v5 = e5->idx, &w5 = e5->w;

            if (!(Data.Reach[v5] & 2) || v1 == v5 || v2 == v5 || v3 == v5) {
              continue;
            } else if (!judge(w4, w5)) {
              continue;
            } else if (v5 == st) {
              if (judge(w5, w1)) {
                ret[2].insert(ret[2].end(), {st, v1, v2, v3, v4});
              }
              continue;
            }

            const DFSEdge *e6 = &G[Head[v5]];

            for (uint it6 = Head[v5]; it6 < Head[v5 + 1]; ++it6, ++e6) {
              const uint &v6 = e6->idx, &w6 = e6->w;

              if (!(Data.Reach[v6] & 1) || v1 == v6 || v2 == v6 || v3 == v6 ||
                  v4 == v6) {
                continue;
              } else if (!judge(w5, w6)) {
                continue;
              } else if (v6 == st) {
                if (judge(w6, w1)) {
                  ret[3].insert(ret[3].end(), {st, v1, v2, v3, v4, v5});
                }
                continue;
              }
              const uint &w7 = Data.LastWeight[v6];

              if (!judge(w6, w7) || !judge(w7, w1)) {
                continue;
              }

              ret[4].insert(ret[4].end(), {st, v1, v2, v3, v4, v5, v6});
            }
          }
        }
      }
    }
  }
}

void HandleFindCycle(uint pid) {
#ifdef TESTSPEED
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif

  auto &Data = ThreadData[pid];

  while (true) {
    while (_JOB_LOCK_.test_and_set())
      ;
    uint job = JobCur < JobsCount ? Jobs[JobCur++] : -1;
    _JOB_LOCK_.clear();
    if (job == -1) break;
    BackSearch(Data, job);
    ForwardSearch(Data, job);
  }

#ifdef TESTSPEED
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ thread %d: [find cycle] [cost: %.4fs]\n", pid, t4 - t1);
#endif
}

void FindCircle() {
  Cycles.reserve(MaxID);
  JobCur = 0;
  std::thread Th[NTHREAD];
  for (uint i = 1; i < NTHREAD; ++i) Th[i] = std::thread(HandleFindCycle, i);
  HandleFindCycle(0);
  for (uint i = 1; i < NTHREAD; ++i) Th[i].join();
}

void CalOffset() {
#ifdef TESTSPEED
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif
  uint tolSize = 0;
  for (uint i = 0; i < JobsCount; ++i) {
    const auto &ret = Cycles[Jobs[i]];
    for (int j = 0; j < 5; ++j) {
      uint sz = ret.cycle[j].size();
      Answers += sz / (j + 3);
      tolSize += sz;
    }
  }

  uint block = tolSize / BUFFERBLOCK + 1;
  uint strow = 0, endrow = 0, stcol = 0, edcol = 0, times = 0;
  for (uint i = 0; i < 5; ++i) {
    for (uint j = 0; j < JobsCount; ++j) {
      if (i == 4 && j == JobsCount - 1) {
        times += Cycles[Jobs[j]].cycle[i].size();
        break;
      }

      if (times >= block) {
        AnswerBuffer[OffSet.size()] = new char[times * 7 * NUMLENGTH];
        OffSet.emplace_back(std::array<uint, 4>{strow, i, stcol, j});
        times = 0;
        strow = i, stcol = j;
      }

      times += Cycles[Jobs[j]].cycle[i].size();
    }
  }

  AnswerBuffer[OffSet.size()] = new char[times * 7 * NUMLENGTH];
  OffSet.emplace_back(std::array<uint, 4>{strow, 4, stcol, JobsCount});
#ifdef TESTSPEED
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ CalOffSet: [offset: %lu] [cost: %.4fs]\n", OffSet.size(), t4 - t1);
#endif
}

void MemCopy(char *&result, const uint &p, const std::vector<uint> &cycle) {
  uint idx = 0;
  for (const auto &v : cycle) {
    const auto &mpid = MapID[v];
    memcpy(result, mpid.str, mpid.len);
    if (++idx == p + 3) {
      idx = 0;
      *(result + mpid.len - 1) = '\n';
    }
    result += mpid.len;
  }
};

void HandleSaveAnswer(uint pid) {
#ifdef TESTSPEED
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif

  while (true) {
    while (_OFFSET_LOCK_.test_and_set())
      ;
    uint cur = OffSetCur < OffSet.size() ? OffSetCur++ : -1;
    _OFFSET_LOCK_.clear();
    if (cur == -1) break;
    const auto &offset = OffSet[cur];
    uint stl = offset[0], edl = offset[1];
    uint stidx = offset[2], edidx = offset[3];
    char *result = AnswerBuffer[cur];
    for (uint i = stl; i <= edl; ++i) {
      uint low = (i == stl ? stidx : 0);
      uint high = (i == edl ? edidx : JobsCount);
      for (uint j = low; j < high; ++j) {
        uint idx = 0;
        const auto &cycle = Cycles[Jobs[j]].cycle[i];
        for (const auto &v : cycle) {
          const auto &mpid = MapID[v];
          memcpy(result, mpid.str, mpid.len);
          if (++idx == i + 3) {
            idx = 0;
            *(result + mpid.len - 1) = '\n';
          }
          result += mpid.len;
        }
      }
    }
    AnswerBufferLen[cur] = result - AnswerBuffer[cur];
    Done[cur] = true;
  }

#ifdef TESTSPEED
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ thread %d: [memcpy] [cost: %.4fs]\n", pid, t4 - t1);
#endif
}

void WriteAnswer(int pid) {
#ifdef TESTSPEED
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif

  sprintf(FirstBuf, "%d\n", Answers);
  FirstBufLen = strlen(FirstBuf);
  FILE *fp = fopen(RESULT, "w");
  fwrite(FirstBuf, 1, FirstBufLen, fp);

  uint idx = 0;
  while (idx < OffSet.size()) {
    while (!Done[idx]) usleep(1);
    fwrite(AnswerBuffer[idx], AnswerBufferLen[idx], 1, fp);
    ++idx;
  }
  fclose(fp);

#ifdef TESTSPEED
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ thread %d: [fwrite] [cost: %.4fs]\n", pid, t4 - t1);
#endif
}

void SaveAnswer() {
  std::thread Th[NTHREAD];
  for (uint i = 1; i < NTHREAD; ++i) {
    Th[i] = std::thread(HandleSaveAnswer, i);
  }
  WriteAnswer(0);
  for (uint i = 1; i < NTHREAD; ++i) Th[i].join();
}

int main() {
  LoadData();
  FindCircle();
  CalOffset();
  SaveAnswer();
  // sleep(1);
#ifdef LOCAL
  std::cerr << "@ Answers: " << Answers << "\n";
#endif
  return 0;
}