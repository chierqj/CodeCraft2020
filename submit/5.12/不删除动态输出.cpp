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
#define TRAIN "../data/19630345/test_data.txt"
#define RESULT "../data/19630345/result.txt"
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

struct DFSEdge {
  uint v, w;
};
struct Edge {
  uint u, v, w;
};
struct PreBuffer {
  uint len;
  char str[NUMLENGTH];
};
struct Answer {
  uint bufsize[5];
  std::vector<uint> cycle[5];
};
struct ThData {
  uint answers = 0;       // 环数目
  uint bufsize = 0;       // bufsize
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
std::atomic_flag _JOB_LOCK_ = ATOMIC_FLAG_INIT;     // job lock
std::atomic_flag _OFFSET_LOCK_ = ATOMIC_FLAG_INIT;  // offset lock
char FirstBuf[NUMLENGTH];                           // 环个数buf
uint ThEdgesCount[NTHREAD];                         // thread: 正向边cnt
uint ThBackEdgesCount[NTHREAD];                     // thread: 反向边cnt
std::vector<std::array<uint, 5>> OffSet;  // 分块输出，每一块位置
uint Jobs[MAXN];                          // 有效点
PreBuffer MapID[MAXN];                    // 解析int
DFSEdge *Children[MAXN][2];               // 子结点
DFSEdge *Parents[MAXN][2];                // 父亲结点
std::vector<Pair> ThChildren[NTHREAD];    // thread: 子结点
std::vector<Pair> ThParents[NTHREAD];     // thread: 父结点
std::vector<Answer> Cycles;               // 所有答案
DFSEdge DFSEdges[MAXEDGE];                // 搜索正向边集合
DFSEdge BackDFSEdges[MAXEDGE];            // 搜索反向边集合
Edge Edges[MAXEDGE];                      // 读入正向边集合
Edge BackEdges[MAXEDGE];                  // 读入反向边集合
Edge ThEdges[NTHREAD][MAXEDGE];           // thread: 正向边集合
Edge ThBackEdges[NTHREAD][MAXEDGE];       // thread: 反向边集合
ThData ThreadData[NTHREAD];               // 线程找环

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
  std::sort(ThEdges[pid], ThEdges[pid] + ThEdgesCount[pid],
            [&](const Edge &e1, const Edge &e2) {
              if (e1.u == e2.u) return e1.v < e2.v;
              return e1.u < e2.u;
            });
  std::sort(ThBackEdges[pid], ThBackEdges[pid] + ThBackEdgesCount[pid],
            [&](const Edge &e1, const Edge &e2) {
              if (e1.v == e2.v) return e1.u > e2.u;
              return e1.v < e2.v;
            });
}

void SortEdge() {
#ifdef LOCAL
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif

  std::thread Th[NTHREAD];
  for (int i = 0; i < NTHREAD; ++i) Th[i] = std::thread(HandleSortEdge, i);
  for (int i = 0; i < NTHREAD; ++i) Th[i].join();
  Edge *ptr = Edges, *ptr1 = BackEdges;
  for (int i = 0; i < NTHREAD; ++i) {
    const uint &cnt = ThEdgesCount[i];
    memcpy(ptr, ThEdges[i], cnt * sizeof(Edge));
    ptr += cnt;
    EdgesCount += cnt;

    const uint &cnt1 = ThBackEdgesCount[i];
    memcpy(ptr1, ThBackEdges[i], cnt1 * sizeof(Edge));
    ptr1 += cnt1;
  }

#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ SortEdge:\t[cost: %.4fs]\n", t4 - t1);
#endif
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
    if (i == 0 || (i > 0 && e.u != pre)) {
      HashID.Insert(e.u);
    }
    HashID.Insert(e.v);
    pre = e.u;
  }
  HashID.Sort();

#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ HashID:\t[cost: %.4fs]\n", t4 - t1);
#endif
}

void HandleCreateSubGraph(int pid) {
  auto &children = ThChildren[pid];
  auto &parents = ThParents[pid];
  parents.reserve(MaxID);
  children.reserve(MaxID);
  for (int i = pid; i < EdgesCount; i += NTHREAD) {
    auto &e = Edges[i];
    const uint &p1 = HashID.Query(e.u);
    const uint &p2 = HashID.Query(e.v);
    DFSEdges[i].v = p2;
    DFSEdges[i].w = e.w;
    if (children[p1].second == 0) {
      children[p1].first = i;
    }
    ++children[p1].second;

    auto &e1 = BackEdges[i];
    const uint &p11 = HashID.Query(e1.u);
    const uint &p22 = HashID.Query(e1.v);
    BackDFSEdges[i].v = p11;
    BackDFSEdges[i].w = e1.w;
    if (parents[p22].second == 0) {
      parents[p22].first = i;
    }
    ++parents[p22].second;
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
    uint minx1 = EdgesCount, len1 = 0;
    for (uint j = 0; j < NTHREAD; ++j) {
      const auto &cdr = ThChildren[j][i];
      if (cdr.second > 0) {
        minx = std::min(minx, cdr.first);
        len += cdr.second;
      }
      const auto &pat = ThParents[j][i];
      if (pat.second > 0) {
        minx1 = std::min(minx1, pat.first);
        len1 += pat.second;
      }
    }
    Children[i][0] = &DFSEdges[minx];
    Children[i][1] = &DFSEdges[minx + len];
    Parents[i][0] = &BackDFSEdges[minx1];
    Parents[i][1] = &BackDFSEdges[minx1 + len1];
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

void addEdge(const uint &u, const uint &v, const uint &w) {
  uint pid = u % NTHREAD;
  auto &e = ThEdges[pid][ThEdgesCount[pid]++];
  e.u = u;
  e.v = v;
  e.w = w;
  pid = v % NTHREAD;
  auto &e1 = ThBackEdges[pid][ThBackEdgesCount[pid]++];
  e1.u = u;
  e1.v = v;
  e1.w = w;
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

  SortEdge();
  CreateHashTable();
  CreateSubGraph();
  CreateGraph();

  for (uint i = 0; i < MaxID; ++i) {
    if (Children[i][1] - Children[i][0] > 0 &&
        Parents[i][1] - Parents[i][0] > 0) {
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
  const auto &f1 = Parents[st];
  for (const auto *e1 = f1[0]; e1 < f1[1]; ++e1) {
    const uint &v1 = e1->v;
    if (v1 <= st) break;
    Data.LastWeight[v1] = e1->w;
    Data.Reach[v1] = 7;
    Data.ReachPoint[Data.ReachCount++] = v1;
    const auto &f2 = Parents[v1];
    for (const auto *e2 = f2[0]; e2 < f2[1]; ++e2) {
      const uint &v2 = e2->v;
      if (v2 <= st) break;
      if (!judge(e2, e1)) continue;
      Data.Reach[v2] |= 6;
      Data.ReachPoint[Data.ReachCount++] = v2;
      const auto &f3 = Parents[v2];
      for (const auto *e3 = f3[0]; e3 < f3[1]; ++e3) {
        const uint &v3 = e3->v;
        if (v3 <= st) break;
        if (v3 == v1 || !judge(e3, e2)) continue;
        Data.Reach[v3] |= 4;
        Data.ReachPoint[Data.ReachCount++] = v3;
      }
    }
  }
}

void ForwardSearch(ThData &Data, const uint &st) {
  uint sz0 = 0, sz1 = 0, sz2 = 0, sz3 = 0, sz4 = 0;

  auto &ret = Cycles[st];
  const auto &c1 = Children[st];
  const uint &len0 = MapID[st].len;

  for (const auto *e1 = c1[0]; e1 < c1[1]; ++e1) {
    if (e1->v < st) continue;

    const uint &len1 = MapID[e1->v].len + len0;
    const auto &c2 = Children[e1->v];

    for (const auto *e2 = c2[0]; e2 < c2[1]; ++e2) {
      if (e2->v <= st || !judge(e1, e2)) continue;

      const uint &len = MapID[e2->v].len + len1;
      const auto &c3 = Children[e2->v];

      for (const auto *e3 = c3[0]; e3 < c3[1]; ++e3) {
        if (e3->v < st || e3->v == e1->v || !judge(e2, e3)) {
          continue;
        } else if (e3->v == st) {
          if (judge(e3, e1)) {
            ret.cycle[0].insert(ret.cycle[0].end(), {st, e1->v, e2->v});
            sz0 += len;
          }
          continue;
        }

        const uint &len3 = MapID[e3->v].len;
        const auto &c4 = Children[e3->v];

        for (const auto *e4 = c4[0]; e4 < c4[1]; ++e4) {
          if (!(Data.Reach[e4->v] & 4) || e1->v == e4->v || e2->v == e4->v) {
            continue;
          } else if (!judge(e3, e4)) {
            continue;
          } else if (e4->v == st) {
            if (judge(e4, e1)) {
              ret.cycle[1].insert(ret.cycle[1].end(),
                                  {st, e1->v, e2->v, e3->v});
              sz1 += len + len3;
            }
            continue;
          }

          const uint &len4 = MapID[e4->v].len;
          const auto &c5 = Children[e4->v];

          for (const auto *e5 = c5[0]; e5 < c5[1]; ++e5) {
            if (!(Data.Reach[e5->v] & 2) || e1->v == e5->v || e2->v == e5->v ||
                e3->v == e5->v) {
              continue;
            } else if (!judge(e4, e5)) {
              continue;
            } else if (e5->v == st) {
              if (judge(e5, e1)) {
                ret.cycle[2].insert(ret.cycle[2].end(),
                                    {st, e1->v, e2->v, e3->v, e4->v});
                sz2 += len + len3 + len4;
              }
              continue;
            }

            const uint &len5 = MapID[e5->v].len;
            const auto &c6 = Children[e5->v];

            for (const auto *e6 = c6[0]; e6 < c6[1]; ++e6) {
              if (!(Data.Reach[e6->v] & 1) || e1->v == e6->v ||
                  e2->v == e6->v || e3->v == e6->v || e4->v == e6->v) {
                continue;
              } else if (!judge(e5, e6)) {
                continue;
              } else if (e6->v == st) {
                if (judge(e6, e1)) {
                  ret.cycle[3].insert(ret.cycle[3].end(),
                                      {st, e1->v, e2->v, e3->v, e4->v, e5->v});
                  sz3 += len + len3 + len4 + len5;
                }
                continue;
              }
              const uint &w7 = Data.LastWeight[e6->v];
              if (!judge(e6->w, w7) || !judge(w7, e1->w)) {
                continue;
              }

              const uint &len6 = MapID[e6->v].len;
              ret.cycle[4].insert(ret.cycle[4].end(), {st, e1->v, e2->v, e3->v,
                                                       e4->v, e5->v, e6->v});
              sz4 += len + len3 + len4 + len5 + len6;
            }
          }
        }
      }
    }
  }
  for (int i = 0; i < 5; ++i) {
    Data.answers += ret.cycle[i].size() / (i + 3);
  }
  Data.bufsize += sz0 + sz1 + sz2 + sz3 + sz4;
  ret.bufsize[0] = sz0;
  ret.bufsize[1] = sz1;
  ret.bufsize[2] = sz2;
  ret.bufsize[3] = sz3;
  ret.bufsize[4] = sz4;
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
  printf("@ thread %d: [find: %d] [cost: %.4fs]\n", pid, Data.answers, t4 - t1);
#endif
}

void FindCircle() {
  Cycles.reserve(MaxID);
  JobCur = 0;
  std::thread Th[NTHREAD];
  for (uint i = 1; i < NTHREAD; ++i) Th[i] = std::thread(HandleFindCycle, i);
  HandleFindCycle(0);
  for (uint i = 1; i < NTHREAD; ++i) Th[i].join();
  for (auto &it : ThreadData) {
    Answers += it.answers;
    TotalBufferSize += it.bufsize;
  }
}

void CalOffset() {
#ifdef TESTSPEED
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif
  uint block = std::sqrt(TotalBufferSize);
  uint strow = 0, endrow = 0, stcol = 0, edcol = 0, times = 0, offsz = 0;
  for (uint i = 0; i < 5; ++i) {
    for (uint j = 0; j < JobsCount; ++j) {
      const uint &job = Jobs[j];
      if (times >= block && !(i == 4 && j == JobsCount - 1)) {
        OffSet.emplace_back(std::array<uint, 5>{strow, i, stcol, j, offsz});
        offsz += times;
        times = 0;
        strow = i, stcol = j;
      }
      times += Cycles[job].bufsize[i];
    }
  }
  OffSet.emplace_back(std::array<uint, 5>{strow, 4, stcol, JobsCount, offsz});
#ifdef TESTSPEED
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ CalOffSet: [offset: %d] [cost: %.4fs]\n", OffSet.size(), t4 - t1);
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

void GetNextOffSet(uint &ret) {
  while (_OFFSET_LOCK_.test_and_set())
    ;
  ret = OffSetCur < OffSet.size() ? OffSetCur++ : -1;
  _OFFSET_LOCK_.clear();
}
void HandleSaveAnswer(uint pid, char *result) {
#ifdef TESTSPEED
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif

  uint offcur = 0;
  char *ptr = result;
  while (true) {
    GetNextOffSet(offcur);
    if (offcur == -1) break;
    const auto &offset = OffSet[offcur];
    uint stl = offset[0], edl = offset[1];
    uint stidx = offset[2], edidx = offset[3];
    uint offsz = FirstBufLen + offset[4];
    result = ptr + offsz;
    for (uint i = stl; i <= edl; ++i) {
      uint low = (i == stl ? stidx : 0);
      uint high = (i == edl ? edidx : JobsCount);
      for (uint j = low; j < high; ++j) {
        WriteAnswer(result, i, Cycles[Jobs[j]].cycle[i]);
      }
    }
  }

#ifdef TESTSPEED
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ thread %d: [save: %.4fs]\n", pid, t4 - t1);
#endif
}

void SaveAnswer() {
  CalOffset();

  sprintf(FirstBuf, "%d\n", Answers);
  FirstBufLen = strlen(FirstBuf);
  TotalBufferSize += FirstBufLen;
  uint fd = open(RESULT, O_RDWR | O_CREAT, 0666);
  char *result = (char *)mmap(NULL, TotalBufferSize, PROT_READ | PROT_WRITE,
                              MAP_SHARED, fd, 0);
  ftruncate(fd, TotalBufferSize);
  close(fd);
  memcpy(result, FirstBuf, FirstBufLen);

  std::thread Th[NTHREAD];
  for (uint i = 1; i < NTHREAD; ++i) {
    Th[i] = std::thread(HandleSaveAnswer, i, result);
  }
  HandleSaveAnswer(0, result);
  for (uint i = 1; i < NTHREAD; ++i) Th[i].join();
}

int main() {
  LoadData();
  FindCircle();
  SaveAnswer();
  sleep(1);
#ifdef LOCAL
  std::cerr << "@ Answers: " << Answers << "";
  std::cerr << ", Bufsize: " << TotalBufferSize << "\n";
#endif
  return 0;
}