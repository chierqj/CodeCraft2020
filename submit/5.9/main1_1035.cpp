#include <fcntl.h>
#include <sys/mman.h>
#include <sys/time.h>
#include <unistd.h>

#include <algorithm>
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
#define TRAIN "./data/19630345/test_data.txt"
#define RESULT "./data/19630345/result.txt"
#else
#define TRAIN "/data/test_data.txt"
#define RESULT "/projects/student/result.txt"
#endif
typedef std::pair<uint, uint> Pair;

/*
 * 常量定义
 */
const uint MAXEDGE = 2000000 + 7;  // 最多边数目
const uint MAXN = MAXEDGE;         // 最多点数目
const uint NTHREAD = 4;            // 线程个数
const uint NUMLENGTH = 12;         // ID最大长度

struct DFSEdge {
  uint v, w;
};
struct Edge {
  uint u, v, w;
  bool operator<(const Edge &r) const {
    if (u == r.u) return v < r.v;
    return u < r.u;
  }
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
} ThreadData[NTHREAD];    // 线程找环

uint MaxID = 0;                                     // 最大点
uint Answers = 0;                                   // 环个数
uint EdgesCount = 0;                                // 边数目
uint JobsCount = 0;                                 // 有效点数目
uint JobCur = 0;                                    // Job光标
std::atomic_flag _JOB_LOCK_ = ATOMIC_FLAG_INIT;     // job lock
uint ThEdgesCount[NTHREAD];                         // count
uint AnswerBufferSize[NTHREAD];                     //
uint OffSet[NTHREAD][4];                            // 偏移量
uint Jobs[MAXN];                                    // 有效点
PreBuffer MapID[MAXN];                              // 解析int
uint Children[MAXN][2];                             // sons
std::vector<Pair> Parents[MAXN];                    // fathers
std::vector<Pair> ThChildren[NTHREAD];              // thsons
std::vector<std::vector<Pair>> ThParents[NTHREAD];  // thfather
std::vector<Answer> Cycles;                         // 所有答案
DFSEdge DFSEdges[MAXEDGE];                          // 所有边
Edge Edges[MAXEDGE];                                // 所有边
Edge ThEdges[NTHREAD][MAXEDGE];                     // 线程边
char *AnswerBuffer[NTHREAD];                        // AnswerBuffer

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
HashTable HashID;

void HandleSortEdge(int pid) {
  std::sort(ThEdges[pid], ThEdges[pid] + ThEdgesCount[pid]);
}
void SortEdgeAndHash() {
#ifdef LOCAL
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif

  std::thread Th[NTHREAD];
  for (int i = 0; i < NTHREAD; ++i) Th[i] = std::thread(HandleSortEdge, i);
  for (int i = 0; i < NTHREAD; ++i) Th[i].join();
  Edge *ptr = Edges;
  for (int i = 0; i < NTHREAD; ++i) {
    const uint &cnt = ThEdgesCount[i];
    memcpy(ptr, ThEdges[i], cnt * sizeof(Edge));
    ptr += cnt;
    EdgesCount += cnt;
  }

#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ SortEdge:\t[cost: %.4fs]\n", t4 - t1);
#endif

#ifdef LOCAL
  struct timeval tim1 {};
  gettimeofday(&tim1, nullptr);
  double t11 = tim1.tv_sec + (tim1.tv_usec / 1000000.0);
#endif

  uint pre = 0;
  for (uint i = 0; i < EdgesCount; ++i) {
    const auto &e = Edges[i];
    if (i == 0 || (i > 0 && e.u != pre)) {
      HashID.Insert(e.u);
    }
    HashID.Insert(e.v);
    pre = e.u;
  }
  HashID.Sort();

#ifdef LOCAL
  gettimeofday(&tim1, nullptr);
  double t41 = tim1.tv_sec + (tim1.tv_usec / 1000000.0);
  printf("@ HashID:\t[cost: %.4fs]\n", t41 - t11);
#endif
}

void HandleCreateSubGraph(int pid) {
  auto &children = ThChildren[pid];
  auto &parents = ThParents[pid];
  parents.reserve(MaxID);
  children.reserve(MaxID);
  for (int i = pid; i < EdgesCount; i += NTHREAD) {
    auto &e = Edges[i];
    e.u = HashID.Query(e.u);
    e.v = HashID.Query(e.v);
    DFSEdges[i].v = e.v;
    DFSEdges[i].w = e.w;
    if (children[e.u].second == 0) {
      children[e.u].first = i;
    }
    ++children[e.u].second;
    parents[e.v].emplace_back(std::make_pair(e.u, e.w));
  }
}
void CreateSubGraph() {
#ifdef LOCAL
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif
  std::thread Th[NTHREAD];
  for (int i = 0; i < NTHREAD; ++i) {
    Th[i] = std::thread(HandleCreateSubGraph, i);
  }
  for (int i = 0; i < NTHREAD; ++i) Th[i].join();
#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ SubGraph:\t[cost: %.4fs]\n", t4 - t1);
#endif
}
void HandleCreateGraph(int pid) {
  for (int i = pid; i < MaxID; i += NTHREAD) {
    uint minx = EdgesCount, len = 0;
    for (uint j = 0; j < NTHREAD; ++j) {
      Parents[i].insert(Parents[i].end(), ThParents[j][i].begin(),
                        ThParents[j][i].end());
      const auto &cdr = ThChildren[j][i];
      if (cdr.second > 0) {
        minx = std::min(minx, cdr.first);
        len += cdr.second;
      }
    }
    // Children[i][0] = &DFSEdges[minx];
    // Children[i][1] = &DFSEdges[minx + len];
    Children[i][0] = minx;
    Children[i][1] = minx + len;
    std::sort(
        Parents[i].begin(), Parents[i].end(),
        [](const Pair &e1, const Pair &e2) { return e1.first > e2.first; });
  }
}
void CreateGraph() {
#ifdef LOCAL
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif
  std::thread Th[NTHREAD];
  for (int i = 0; i < NTHREAD; ++i) Th[i] = std::thread(HandleCreateGraph, i);
  for (int i = 0; i < NTHREAD; ++i) Th[i].join();
#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ Graph:\t[cost: %.4fs]\n", t4 - t1);
#endif
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
    uint pid = u % NTHREAD;
    auto &cnt = ThEdgesCount[pid];
    auto &e = ThEdges[pid][cnt++];
    e.u = u;
    e.v = v;
    e.w = w;
    u = v = w = 0;
  }

#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ Buffer:\t[cost: %.4fs]\n", t4 - t1);
#endif

  SortEdgeAndHash();
  CreateSubGraph();
  CreateGraph();

  for (uint i = 0; i < MaxID; ++i) {
    if (Children[i][1] > 0 && !Parents[i].empty()) {
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

inline bool judge(const DFSEdge *e1, const DFSEdge *e2) {
  return (e2->w > W5_MAX || e1->w <= P5(e2->w)) &&
         (e1->w > W3_MAX || P3(e1->w) >= e2->w);
}

void BackSearch(ThData &Data, const uint &st) {
  for (uint i = 0; i < Data.ReachCount; ++i) {
    const uint &v = Data.ReachPoint[i];
    Data.Reach[v] = 0;
  }
  Data.ReachCount = 0;
  Data.ReachPoint[Data.ReachCount++] = st;
  Data.Reach[st] = 7;
  const auto &parent1 = Parents[st];
  for (const auto &it1 : parent1) {
    const uint &v1 = it1.first;
    if (v1 <= st) break;
    const uint &w1 = it1.second;
    Data.LastWeight[v1] = w1;
    Data.Reach[v1] = 7;
    Data.ReachPoint[Data.ReachCount++] = v1;
    const auto &parent2 = Parents[v1];
    for (const auto &it2 : parent2) {
      const uint &v2 = it2.first;
      if (v2 <= st) break;
      const uint &w2 = it2.second;
      if (!judge(w2, w1)) continue;
      Data.Reach[v2] |= 6;
      Data.ReachPoint[Data.ReachCount++] = v2;
      const auto &parent3 = Parents[v2];
      for (const auto &it3 : parent3) {
        const uint &v3 = it3.first;
        if (v3 <= st) break;
        const uint &w3 = it3.second;
        if (v3 == v1) continue;
        if (!judge(w3, w2)) continue;
        Data.Reach[v3] |= 4;
        Data.ReachPoint[Data.ReachCount++] = v3;
      }
    }
  }
}

/*
void ForwardSearch(ThData &Data, const uint &st) {
  auto &ret = Cycles[st];
  const auto &c1 = Children[st];
  for (const DFSEdge *e1 = c1[0]; e1 < c1[1]; ++e1) {
    const uint &v1 = e1->v;
    if (v1 < st) continue;
    const uint &w1 = e1->w;
    const auto &c2 = Children[v1];
    for (const DFSEdge *e2 = c2[0]; e2 < c2[1]; ++e2) {
      const uint &v2 = e2->v;
      if (v2 <= st || !judge(e1, e2)) continue;
      const auto &c3 = Children[v2];
      for (const DFSEdge *e3 = c3[0]; e3 < c3[1]; ++e3) {
        const uint &v3 = e3->v;
        if (v3 < st || v3 == v1 || !judge(e2, e3)) {
          continue;
        } else if (v3 == st) {
          if (!judge(e3, e1)) continue;
          ret.cycle[0].insert(ret.cycle[0].end(), {st, v1, v2});
          continue;
        }
        const auto &c4 = Children[v3];
        for (const DFSEdge *e4 = c4[0]; e4 < c4[1]; ++e4) {
          const uint &v4 = e4->v;
          if (!(Data.Reach[v4] & 4) || v1 == v4 || v2 == v4) {
            continue;
          } else if (!judge(e3, e4)) {
            continue;
          } else if (v4 == st) {
            if (!judge(e4, e1)) continue;
            ret.cycle[1].insert(ret.cycle[1].end(), {st, v1, v2, v3});
            continue;
          }
          const auto &c5 = Children[v4];
          for (const DFSEdge *e5 = c5[0]; e5 < c5[1]; ++e5) {
            const uint &v5 = e5->v;
            if (!(Data.Reach[v5] & 2) || v1 == v5 || v2 == v5 || v3 == v5) {
              continue;
            } else if (!judge(e4, e5)) {
              continue;
            } else if (v5 == st) {
              if (!judge(e5, e1)) continue;
              ret.cycle[2].insert(ret.cycle[2].end(), {st, v1, v2, v3, v4});
              continue;
            }
            const auto &c6 = Children[v5];
            for (const DFSEdge *e6 = c6[0]; e6 < c6[1]; ++e6) {
              const uint &v6 = e6->v;
              if (!(Data.Reach[v6] & 1) || v1 == v6 || v2 == v6 || v3 == v6 ||
                  v4 == v6) {
                continue;
              } else if (!judge(e5, e6)) {
                continue;
              } else if (v6 == st) {
                if (!judge(e6, e1)) continue;
                ret.cycle[3].insert(ret.cycle[3].end(),
                                    {st, v1, v2, v3, v4, v5});
                continue;
              }
              const uint &w7 = Data.LastWeight[v6];
              if (!judge(e6->w, w7) || !judge(w7, w1)) {
                continue;
              }
              ret.cycle[4].insert(ret.cycle[4].end(),
                                  {st, v1, v2, v3, v4, v5, v6});
            }
          }
        }
      }
    }
  }
}
*/
void ForwardSearch(ThData &Data, const uint &st) {
  auto &ret = Cycles[st];
  const auto &cdr1 = Children[st];
  const DFSEdge *e1 = &DFSEdges[cdr1[0]];
  for (uint it1 = cdr1[0]; it1 != cdr1[1]; ++it1, ++e1) {
    const uint &v1 = e1->v;
    if (v1 < st) continue;
    const auto &cdr2 = Children[v1];
    const DFSEdge *e2 = &DFSEdges[cdr2[0]];
    for (uint it2 = cdr2[0]; it2 != cdr2[1]; ++it2, ++e2) {
      const uint &v2 = e2->v;
      if (v2 <= st || !judge(e1, e2)) continue;
      const auto &cdr3 = Children[v2];
      const DFSEdge *e3 = &DFSEdges[cdr3[0]];
      for (uint it3 = cdr3[0]; it3 != cdr3[1]; ++it3, ++e3) {
        const uint &v3 = e3->v;
        if (v3 < st || v3 == v1 || !judge(e2, e3)) {
          continue;
        } else if (v3 == st) {
          if (!judge(e3, e1)) continue;
          ret.cycle[0].insert(ret.cycle[0].end(), {st, v1, v2});
          continue;
        }
        const auto &cdr4 = Children[v3];
        const DFSEdge *e4 = &DFSEdges[cdr4[0]];
        for (uint it4 = cdr4[0]; it4 != cdr4[1]; ++it4, ++e4) {
          const uint &v4 = e4->v;
          if (!(Data.Reach[v4] & 4) || v1 == v4 || v2 == v4) {
            continue;
          } else if (!judge(e3, e4)) {
            continue;
          } else if (v4 == st) {
            if (!judge(e4, e1)) continue;
            ret.cycle[1].insert(ret.cycle[1].end(), {st, v1, v2, v3});
            continue;
          }
          const auto &cdr5 = Children[v4];
          const DFSEdge *e5 = &DFSEdges[cdr5[0]];
          for (uint it5 = cdr5[0]; it5 != cdr5[1]; ++it5, ++e5) {
            const uint &v5 = e5->v;
            if (!(Data.Reach[v5] & 2) || v1 == v5 || v2 == v5 || v3 == v5) {
              continue;
            } else if (!judge(e4, e5)) {
              continue;
            } else if (v5 == st) {
              if (!judge(e5, e1)) continue;
              ret.cycle[2].insert(ret.cycle[2].end(), {st, v1, v2, v3, v4});
              continue;
            }
            const auto &cdr6 = Children[v5];
            const DFSEdge *e6 = &DFSEdges[cdr6[0]];
            for (uint it6 = cdr6[0]; it6 != cdr6[1]; ++it6, ++e6) {
              const uint &v6 = e6->v;
              if (!(Data.Reach[v6] & 1) || v1 == v6 || v2 == v6 || v3 == v6 ||
                  v4 == v6) {
                continue;
              } else if (!judge(e5, e6)) {
                continue;
              } else if (v6 == st) {
                if (!judge(e6, e1)) continue;
                ret.cycle[3].insert(ret.cycle[3].end(),
                                    {st, v1, v2, v3, v4, v5});
                continue;
              }
              const uint &w7 = Data.LastWeight[v6];
              if (!judge(e6->w, w7) || !judge(w7, e1->w)) {
                continue;
              }
              ret.cycle[4].insert(ret.cycle[4].end(),
                                  {st, v1, v2, v3, v4, v5, v6});
            }
          }
        }
      }
    }
  }
}

inline void GetNextJob(uint &job) {
  while (_JOB_LOCK_.test_and_set())
    ;
  job = JobCur < JobsCount ? Jobs[JobCur++] : -1;
  _JOB_LOCK_.clear();
}

void HandleFindCycle(uint pid) {
#ifdef TESTSPEED
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif

  uint job = 0;
  auto &Data = ThreadData[pid];

  while (true) {
    GetNextJob(job);
    if (job == -1) break;
    BackSearch(Data, job);
    ForwardSearch(Data, job);
  }

#ifdef TESTSPEED
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ thread %d: FindCycle [cost: %.4fs]\n", pid, t4 - t1);
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

  uint TolSize = 0;
  for (int i = 0; i < JobsCount; ++i) {
    const uint &job = Jobs[i];
    const auto &ans = Cycles[job];
    for (int j = 0; j < 5; ++j) {
      const uint &sz = ans.cycle[j].size();
      TolSize += sz;
      Answers += sz / (j + 3);
    }
  }
  uint block = TolSize / NTHREAD;
  uint x = 0, stl = 0, stidx = 0;
  uint tidx = 0;
  for (uint i = 0; i < 5; ++i) {
    for (uint j = 0; j < JobsCount; ++j) {
      const uint &sz = Cycles[Jobs[j]].cycle[i].size();
      if (x < block) {
        x += sz;
      } else {
        AnswerBuffer[tidx] = new char[x * NUMLENGTH];
        OffSet[tidx][0] = stl;
        OffSet[tidx][1] = i;
        OffSet[tidx][2] = stidx;
        OffSet[tidx++][3] = j;
        stl = i;
        stidx = j;
        x = 0;
      }
    }
  }
  AnswerBuffer[tidx] = new char[x * NUMLENGTH];
  OffSet[tidx][0] = stl;
  OffSet[tidx][1] = 4;
  OffSet[tidx][2] = stidx;
  OffSet[tidx][3] = JobsCount;

#ifdef TESTSPEED
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ CalOffSet [cost: %.4fs]\n", t4 - t1);
#endif
}

void WriteAnswer(char *&result, const uint &p, const std::vector<uint> &cycle) {
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
  char *result = AnswerBuffer[pid];
  uint job = 0;
  const auto &offset = OffSet[pid];
  uint stl = offset[0], edl = offset[1];
  uint stidx = offset[2], edidx = offset[3];
  for (uint i = stl; i <= edl; ++i) {
    uint low = (i == stl ? stidx : 0);
    uint high = (i == edl ? edidx : JobsCount);
    for (uint j = low; j < high; ++j) {
      uint idx = 0;
      for (const auto &v : Cycles[Jobs[j]].cycle[i]) {
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
  AnswerBufferSize[pid] = result - AnswerBuffer[pid];

#ifdef TESTSPEED
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ thread %d: CreateBuffer [cost: %.4fs]\n", pid, t4 - t1);
#endif
}

void SaveAnswer() {
  std::thread Th[NTHREAD];
  for (uint i = 0; i < NTHREAD; ++i) Th[i] = std::thread(HandleSaveAnswer, i);
  for (uint i = 0; i < NTHREAD; ++i) Th[i].join();

#ifdef TESTSPEED
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif
  char tmp[NUMLENGTH];
  sprintf(tmp, "%d\n", Answers);
  uint firlen = (uint)(std::log10(Answers) + 2);
  FILE *fp = fopen(RESULT, "w");
  fwrite(tmp, 1, firlen, fp);
  for (int i = 0; i < NTHREAD; ++i) {
    fwrite(AnswerBuffer[i], 1, AnswerBufferSize[i], fp);
  }
  fclose(fp);
#ifdef TESTSPEED
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ Fwrite: [cost: %.4fs]\n", t4 - t1);
#endif
}

int main() {
  LoadData();
  FindCircle();
  CalOffset();
  SaveAnswer();
#ifdef LOCAL
  std::cerr << "@ Answers: " << Answers << "\n";
#endif
  return 0;
}