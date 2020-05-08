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
const uint MAXN = MAXEDGE << 1;    // 最多点数目
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
  uint bufsize[5];
  std::vector<uint> cycle0;
  std::vector<uint> cycle1;
  std::vector<uint> cycle2;
  std::vector<uint> cycle3;
  std::vector<uint> cycle4;
};
struct ThData {
  uint answers = 0;       // 环数目
  uint bufsize = 0;       // bufsize
  uint ReachCount = 0;    // 反向可达点数目
  char Reach[MAXN];       // 标记反向可达
  uint ReachPoint[MAXN];  // 可达点集合
  uint LastWeight[MAXN];  // 最后一步权重
} ThreadData[NTHREAD];    // 线程找环

uint MaxID = 0;                                     // 最大点
uint Answers = 0;                                   // 环个数
uint EdgesCount = 0;                                // 边数目
uint JobsCount = 0;                                 // 有效点数目
uint TotalBufferSize = 0;                           // 总buffer大小
uint FirstBufLen = 0;                               // 换个数bufsize
uint JobCur = 0;                                    // Job光标
std::atomic_flag _JOB_LOCK_ = ATOMIC_FLAG_INIT;     // job lock
char FirstBuf[NUMLENGTH];                           // 环个数buf
uint ThEdgesCount[NTHREAD];                         // count
uint OffSet[NTHREAD][5];                            // 偏移量
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
      mpid.len = (uint)(std::log10(Map[idx].key) + 2);
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

void ForwardSearch(ThData &Data, const uint &st) {
  uint sz0 = 0, sz1 = 0, sz2 = 0, sz3 = 0, sz4 = 0;
  const uint &len0 = MapID[st].len;
  auto &ret = Cycles[st];
  const auto &cdr1 = Children[st];
  const DFSEdge *e1 = &DFSEdges[cdr1[0]];
  for (uint it1 = cdr1[0]; it1 != cdr1[1]; ++it1, ++e1) {
    if (e1->v < st) continue;
    const uint &len1 = MapID[e1->v].len + len0;
    const auto &cdr2 = Children[e1->v];
    const DFSEdge *e2 = &DFSEdges[cdr2[0]];
    for (uint it2 = cdr2[0]; it2 != cdr2[1]; ++it2, ++e2) {
      if (e2->v <= st || !judge(e1, e2)) continue;
      const uint &len = MapID[e2->v].len + len1;
      const auto &cdr3 = Children[e2->v];
      const DFSEdge *e3 = &DFSEdges[cdr3[0]];
      for (uint it3 = cdr3[0]; it3 != cdr3[1]; ++it3, ++e3) {
        if (e3->v < st || e3->v == e1->v || !judge(e2, e3)) {
          continue;
        } else if (e3->v == st) {
          if (!judge(e3, e1)) continue;
          ret.cycle0.insert(ret.cycle0.end(), {st, e1->v, e2->v});
          sz0 += len;
          continue;
        }
        const uint &len3 = MapID[e3->v].len;
        const auto &cdr4 = Children[e3->v];
        const DFSEdge *e4 = &DFSEdges[cdr4[0]];
        for (uint it4 = cdr4[0]; it4 != cdr4[1]; ++it4, ++e4) {
          if (!(Data.Reach[e4->v] & 4) || e1->v == e4->v || e2->v == e4->v) {
            continue;
          } else if (!judge(e3, e4)) {
            continue;
          } else if (e4->v == st) {
            if (!judge(e4, e1)) continue;
            ret.cycle1.insert(ret.cycle1.end(), {st, e1->v, e2->v, e3->v});
            sz1 += len + len3;
            continue;
          }
          const uint &len4 = MapID[e4->v].len;
          const auto &cdr5 = Children[e4->v];
          const DFSEdge *e5 = &DFSEdges[cdr5[0]];
          for (uint it5 = cdr5[0]; it5 != cdr5[1]; ++it5, ++e5) {
            if (!(Data.Reach[e5->v] & 2) || e1->v == e5->v || e2->v == e5->v ||
                e3->v == e5->v) {
              continue;
            } else if (!judge(e4, e5)) {
              continue;
            } else if (e5->v == st) {
              if (!judge(e5, e1)) continue;
              ret.cycle2.insert(ret.cycle2.end(),
                                {st, e1->v, e2->v, e3->v, e4->v});
              sz2 += len + len3 + len4;
              continue;
            }
            const uint &len5 = MapID[e5->v].len;
            const auto &cdr6 = Children[e5->v];
            const DFSEdge *e6 = &DFSEdges[cdr6[0]];
            for (uint it6 = cdr6[0]; it6 != cdr6[1]; ++it6, ++e6) {
              if (!(Data.Reach[e6->v] & 1) || e1->v == e6->v ||
                  e2->v == e6->v || e3->v == e6->v || e4->v == e6->v) {
                continue;
              } else if (!judge(e5, e6)) {
                continue;
              } else if (e6->v == st) {
                if (!judge(e6, e1)) continue;
                ret.cycle3.insert(ret.cycle3.end(),
                                  {st, e1->v, e2->v, e3->v, e4->v, e5->v});
                sz3 += len + len3 + len4 + len5;
                continue;
              }
              const uint &w7 = Data.LastWeight[e6->v];
              if (!judge(e6->w, w7) || !judge(w7, e1->w)) {
                continue;
              }
              const uint &len6 = MapID[e6->v].len;
              ret.cycle4.insert(ret.cycle4.end(),
                                {st, e1->v, e2->v, e3->v, e4->v, e5->v, e6->v});
              sz4 += len + len3 + len4 + len5 + len6;
            }
          }
        }
      }
    }
  }
  Data.answers +=
      (ret.cycle0.size() / 3 + ret.cycle1.size() / 4 + ret.cycle2.size() / 5 +
       ret.cycle3.size() / 6 + ret.cycle4.size() / 7);
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
  uint block = TotalBufferSize / NTHREAD;
  uint tol = 0, x = 0, stl = 0, stidx = 0;
  uint tidx = 0;
  for (uint i = 0; i < 5; ++i) {
    for (uint j = 0; j < JobsCount; ++j) {
      if (x > block) {
        OffSet[tidx][0] = stl;
        OffSet[tidx][1] = i;
        OffSet[tidx][2] = stidx;
        OffSet[tidx][3] = j;
        OffSet[tidx++][4] = tol;
        stl = i;
        stidx = j;
        tol += x;
        x = 0;
      }
      x += Cycles[Jobs[j]].bufsize[i];
    }
  }
  OffSet[tidx][0] = stl;
  OffSet[tidx][1] = 4;
  OffSet[tidx][2] = stidx;
  OffSet[tidx][3] = JobsCount;
  OffSet[tidx++][4] = tol;
  // OffSet[tidx++] = std::make_tuple(stl, 4, stidx, JobsCount, tol);
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

void HandleSaveAnswer(uint pid, char *result) {
#ifdef TESTSPEED
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif

  char *ptr = result;
  uint job = 0;
  const auto &offset = OffSet[pid];
  uint stl = offset[0], edl = offset[1];
  uint stidx = offset[2], edidx = offset[3];
  uint offsz = FirstBufLen + offset[4];
  result = ptr + offsz;

  for (uint i = stl; i <= edl; ++i) {
    uint low = (i == stl ? stidx : 0);
    uint high = (i == edl ? edidx : JobsCount);
    switch (i) {
      case 0:
        for (uint j = low; j < high; ++j) {
          WriteAnswer(result, i, Cycles[Jobs[j]].cycle0);
        }
        break;
      case 1:
        for (uint j = low; j < high; ++j) {
          WriteAnswer(result, i, Cycles[Jobs[j]].cycle1);
        }
        break;
      case 2:
        for (uint j = low; j < high; ++j) {
          WriteAnswer(result, i, Cycles[Jobs[j]].cycle2);
        }
        break;
      case 3:
        for (uint j = low; j < high; ++j) {
          WriteAnswer(result, i, Cycles[Jobs[j]].cycle3);
        }
        break;
      case 4:
        for (uint j = low; j < high; ++j) {
          WriteAnswer(result, i, Cycles[Jobs[j]].cycle4);
        }
        break;
      default:
        break;
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

  JobCur = 0;

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
#ifdef LOCAL
  std::cerr << "@ Answers: " << Answers << "";
  std::cerr << ", Bufsize: " << TotalBufferSize << "\n";
#endif
  return 0;
}