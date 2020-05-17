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
#define ulong uint64_t
#define W3_MAX 715827882
#define W5_MAX 429496729
#define INF 4294967295
#define P3(x) ((x << 1) + x)
#define P5(x) ((x << 2) + x)
#define P10(x) ((x << 3) + (x << 1))

#ifdef LOCAL
#define TRAIN "../data/15371869/test_data.txt"
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
const uint MAXN = MAXEDGE;         // 最多点数目
const uint NTHREAD = 4;            // 线程个数
const uint NUMLENGTH = 12;         // ID最大长度
const uint BUFFERBLOCK = 64;       // 输出

struct DFSEdge {
  uint idx;
  ulong w;
};
struct Edge {
  uint u, v;
  ulong w;
};
struct PreBuffer {
  uint len;
  char str[NUMLENGTH];
};
struct Answer {
  std::vector<uint> cycle[6];
};
struct ThData1 {
  uint ReachCount = 0;     // 反向可达点数目
  uint Reach[MAXN];        // 标记反向可达
  uint ReachPoint[MAXN];   // 可达点集合
  ulong LastWeight[MAXN];  // 最后一步权重
};

struct P2 {
  uint p1;
  ulong w;
};
struct P3 {
  uint p1, p2;
  ulong w;
};
struct P4 {
  uint p1, p2, p3;
  ulong w;
};

struct ThData2 {
  uint ReachCount = 0;     // 反向可达点数目
  uint Reach[MAXN];        // 标记反向可达
  bool vis[MAXN];          // 是否sort
  uint ReachPoint[MAXN];   // 可达点集合
  ulong LastWeight[MAXN];  // 最后一步权重
  std::vector<std::vector<P2>> Path2;
  std::vector<std::vector<P3>> Path3;
  std::vector<std::vector<P4>> Path4;
};

uint MaxID = 0;                                     // 最大点
uint Answers = 0;                                   // 环个数
uint EdgesCount = 0;                                // 边数目
uint JobsCount = 0;                                 // 有效点数目
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
std::vector<std::array<uint, 4>> OffSet;            // 分块输出
uint Jobs[MAXN];                                    // 有效点
PreBuffer MapID[MAXN];                              // 解析int
std::vector<Answer> Cycles;                         // 所有答案
ThData1 ThreadData1[NTHREAD];                       // 线程找环
ThData2 ThreadData2[NTHREAD];                       // 线程找环
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

  uint Size() { return count; }
  uint hash(const uint &k, const int &i) {
    return (k % MOD1 + i * (MOD2 - k % MOD2)) % MOD1;
  }
  void Insert(const uint &key) {
    for (int i = 0; i < MOD1; ++i) {
      const uint &val = this->hash(key, i);
      if (Map[val].val == -1) {
        Map[val].key = key;
        Map[val].val = count;
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
  // void Sort() {
  //   std::sort(HashIdx, HashIdx + count, [&](const uint &x, const uint &y) {
  //     return Map[x].key < Map[y].key;
  //   });
  //   for (int i = 0; i < count; ++i) {
  //     const uint &idx = HashIdx[i];
  //     Map[idx].val = i;
  //     auto &mpid = MapID[i];
  //     sprintf(mpid.str, "%d,", Map[idx].key);
  //     mpid.len = strlen(mpid.str);
  //   }
  // }
};
HashTable HashMap;

struct LoadInfo {
  uint edgeCount = 0;
  Edge edges[MAXEDGE];
  uint thedgesCnt[NTHREAD];
  Edge thedges[NTHREAD][MAXEDGE / NTHREAD];
  uint offsz = 0;
  HashTable hashmap;
};
LoadInfo LoadInfos[NTHREAD];

void addEdge(const uint &u, const uint &v, const ulong &w, LoadInfo &data) {
  uint mod = u % NTHREAD;
  auto &e = data.thedges[mod][data.thedgesCnt[mod]];
  e.u = u;
  e.v = v;
  e.w = w;
  ++data.thedgesCnt[mod];
}
void HandleReadBuffer(const char *buffer, uint st, uint ed, uint pid) {
  const char *ptr = buffer + st, *end = buffer + ed;
  uint u = 0, v = 0;
  ulong w = 0;
  auto &loadinfo = LoadInfos[pid];
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

    bool f = false;
    uint cnt = 0;
    while (*ptr != '\r' && *ptr != '\n') {
      if (*ptr == '.') {
        ++ptr;
        f = true;
        continue;
      }
      if (f) {
        ++cnt;
      }
      w = w * 10 + *ptr - '0';
      ++ptr;
    }
    if (*ptr == '\r') ++ptr;
    ++ptr;
    if (cnt == 0) w = w * 100;
    if (cnt == 1) w = w * 10;
    addEdge(u, v, w, loadinfo);
    u = v = w = 0;
  }
}
void ReadBuffer() {
  uint fd = open(TRAIN, O_RDONLY);
  uint bufsize = lseek(fd, 0, SEEK_END);
  char *buffer = (char *)mmap(NULL, bufsize, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);

  std::thread Th[NTHREAD];
  uint st = 0, block = bufsize / NTHREAD;
  for (uint i = 0; i < NTHREAD; ++i) {
    if (i == NTHREAD - 1) {
      Th[i] = std::thread(HandleReadBuffer, buffer, st, bufsize, i);
      break;
    }
    uint ed = st + block;
    while (buffer[ed] != '\n') ++ed;
    ++ed;
    Th[i] = std::thread(HandleReadBuffer, buffer, st, ed, i);
    st = ed;
  }
  for (auto &it : Th) it.join();
}
/*
 * LoadInfos[pid]: edges只存%4的边,，所以需要memcpy
 */
void SortEdge(uint pid) {
  auto &data = LoadInfos[pid];
  Edge *ptr = data.edges;
  for (uint i = 0; i < NTHREAD; ++i) {
    const auto &info = LoadInfos[i];
    memcpy(ptr, info.thedges[pid], info.thedgesCnt[pid] * sizeof(Edge));
    ptr += info.thedgesCnt[pid];
    data.edgeCount += info.thedgesCnt[pid];
  }
  std::sort(data.edges, data.edges + data.edgeCount,
            [&](const Edge &e1, const Edge &e2) {
              if (e1.u == e2.u) return e1.v < e2.v;
              return e1.u < e2.u;
            });
}
void HashEdge(uint pid) {
  auto &data = LoadInfos[pid];
  for (uint i = 0; i < data.edgeCount; ++i) {
    const auto &e = data.edges[i];
    data.hashmap.Insert(e.u);
  }
}
std::vector<uint> IDDom;
void CalMapOffSet() {
  for (uint i = 1; i < NTHREAD; ++i) {
    const uint &sz = LoadInfos[i - 1].hashmap.Size();
    LoadInfos[i].offsz = LoadInfos[i - 1].offsz + sz;
  }
  MaxID = LoadInfos[NTHREAD - 1].offsz + LoadInfos[NTHREAD - 1].hashmap.Size();
  IDDom.reserve(MaxID);
}
void MergeHashTable(uint pid) {
  auto &data = LoadInfos[pid];
  uint offsz = LoadInfos[pid].offsz;
  for (uint i = 0; i < data.hashmap.Size(); ++i) {
    auto &p = data.hashmap.Map[data.hashmap.HashIdx[i]];
    p.val += offsz;
    IDDom[p.val] = p.key;
  }
}
std::vector<uint> Rank;
void RankHashTable() {
  std::vector<uint> vec(MaxID);
  Rank.reserve(MaxID);
  for (uint i = 0; i < MaxID; ++i) vec[i] = i;
  std::sort(vec.begin(), vec.end(),
            [&](const uint &x, const uint &y) { return IDDom[x] < IDDom[y]; });
  for (uint i = 0; i < MaxID; ++i) {
    const uint &x = vec[i];
    Rank[x] = i;
  }
}
void ReHashTable(uint pid) {
  auto &hashmap = LoadInfos[pid].hashmap;
  for (uint i = 0; i < hashmap.Size(); ++i) {
    auto &p = hashmap.Map[hashmap.HashIdx[i]];
    p.val = Rank[p.val];
    auto &mpid = MapID[p.val];
    sprintf(mpid.str, "%d,", p.key);
    mpid.len = strlen(mpid.str);
  }
}
void RankEdge(uint pid) {
  auto &data = LoadInfos[pid];
  for (uint i = 0; i < data.edgeCount; ++i) {
    auto &e = data.edges[i];
    e.u = data.hashmap.Query(e.u);
    int p = LoadInfos[e.v % NTHREAD].hashmap.Query(e.v);
    e.v = (p == -1 ? INF : p);
  }
}

Edge Edges[MAXEDGE];
void MergeEdge() {
  uint left[NTHREAD] = {0};
  Edge *ptr = Edges;
  while (true) {
    uint minx = MaxID + 7, minidx = -1;
    for (uint i = 0; i < NTHREAD; ++i) {
      const auto &data = LoadInfos[i];
      if (left[i] < data.edgeCount && data.edges[left[i]].u < minx) {
        minx = data.edges[left[i]].u;
        minidx = i;
      }
    }
    if (minidx == -1) break;
    const auto &data = LoadInfos[minidx];
    uint &l = left[minidx];
    while (data.edges[l].u == minx) {
      if (data.edges[l].v != INF) {
        memcpy(ptr, &data.edges[left[minidx]], sizeof(Edge));
        ++ptr;
      }
      ++l;
    }
  }
  EdgesCount = ptr - Edges;
}
void BuildGraph(uint pid) {
  uint pre = INF;
  for (uint i = 0; i < EdgesCount; ++i) {
    const auto &e = Edges[i];
    if (e.u % NTHREAD == pid) {
      if (e.u != pre) {
        Head[e.u] = i;
      }
      ++HeadLen[e.u];
      pre = e.u;
      G[i].idx = e.v;
      G[i].w = e.w;
    }
    if (e.v % NTHREAD == pid) {
      ++BackLen[e.v];
    }
  }
}
void BuildBackGraph(uint pid) {
  std::vector<uint> cnt(MaxID, 0);
  for (int i = EdgesCount - 1; i >= 0; --i) {
    const auto &e = Edges[i];
    if (e.v % NTHREAD == pid) {
      auto &p = GBack[Back[e.v] + cnt[e.v]];
      p.idx = e.u;
      p.w = e.w;
      ++cnt[e.v];
    }
  }
}

void LoadData() {
#ifdef LOCAL
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif
  ReadBuffer();

  std::thread Th[NTHREAD];

  for (uint i = 0; i < NTHREAD; ++i) Th[i] = std::thread(SortEdge, i);
  for (auto &it : Th) it.join();

  for (uint i = 0; i < NTHREAD; ++i) Th[i] = std::thread(HashEdge, i);
  for (auto &it : Th) it.join();

  CalMapOffSet();

  for (uint i = 0; i < NTHREAD; ++i) Th[i] = std::thread(MergeHashTable, i);
  for (auto &it : Th) it.join();

  RankHashTable();

  for (uint i = 0; i < NTHREAD; ++i) Th[i] = std::thread(ReHashTable, i);
  for (auto &it : Th) it.join();

  for (uint i = 0; i < NTHREAD; ++i) Th[i] = std::thread(RankEdge, i);
  for (auto &it : Th) it.join();

  MergeEdge();

  for (uint i = 0; i < NTHREAD; ++i) Th[i] = std::thread(RankEdge, i);
  for (auto &it : Th) it.join();

  for (uint i = 0; i < NTHREAD; ++i) Th[i] = std::thread(BuildGraph, i);
  for (auto &it : Th) it.join();

  for (uint i = 1; i <= MaxID; ++i) {
    Head[i] = Head[i - 1] + HeadLen[i - 1];
    Back[i] = Back[i - 1] + BackLen[i - 1];
  }

  for (uint i = 0; i < NTHREAD; ++i) Th[i] = std::thread(BuildBackGraph, i);
  for (auto &it : Th) it.join();

  for (uint i = 0; i < MaxID; ++i) {
    if (HeadLen[i] > 0 && BackLen[i] > 0) {
      Jobs[JobsCount++] = i;
    }
  }

  // for (uint i = 0; i < EdgesCount; ++i) {
  //   auto &e = Edges[i];
  //   std::cerr << e.u << ", " << e.v << ", " << e.w << "\n";
  // }

#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  std::cerr << "MaxID: " << MaxID << ", E: " << EdgesCount << "\n";
  printf("@ LoadData:\t[cost: %.4fs]\n", t4 - t1);
#endif
}

// inline bool judge(const uint &w1, const uint &w2) {
//   return (w2 > W5_MAX || w1 <= P5(w2)) && (w1 > W3_MAX || P3(w1) >= w2);
// }
inline bool judge(const ulong &w1, const ulong &w2) {
  return (w1 <= P5(w2)) && (P3(w1) >= w2);
}
// 8 4 2 1
void BackSearch1(ThData1 &Data, const uint &st) {
  for (uint i = 0; i < Data.ReachCount; ++i) {
    const uint &v = Data.ReachPoint[i];
    Data.Reach[v] = 0;
  }
  Data.ReachCount = 0;
  Data.ReachPoint[Data.ReachCount++] = st;
  Data.Reach[st] = 15;

  const DFSEdge *e1 = &GBack[Back[st]];
  for (uint it1 = Back[st]; it1 < Back[st + 1]; ++it1, ++e1) {
    const uint &v1 = e1->idx;
    if (v1 <= st) break;
    const ulong &w1 = e1->w;
    Data.LastWeight[v1] = w1;
    Data.Reach[v1] = 15;
    Data.ReachPoint[Data.ReachCount++] = v1;

    const DFSEdge *e2 = &GBack[Back[v1]];
    for (uint it2 = Back[v1]; it2 < Back[v1 + 1]; ++it2, ++e2) {
      const uint &v2 = e2->idx;
      if (v2 <= st) break;
      const ulong &w2 = e2->w;

      if (!judge(w2, w1)) continue;
      Data.Reach[v2] |= 14;
      Data.ReachPoint[Data.ReachCount++] = v2;

      const DFSEdge *e3 = &GBack[Back[v2]];
      for (uint it3 = Back[v2]; it3 < Back[v2 + 1]; ++it3, ++e3) {
        const uint &v3 = e3->idx;
        if (v3 <= st) break;
        if (v3 == v1) continue;
        const ulong &w3 = e3->w;
        if (!judge(w3, w2)) continue;
        Data.Reach[v3] |= 12;
        Data.ReachPoint[Data.ReachCount++] = v3;
        const DFSEdge *e4 = &GBack[Back[v3]];
        for (uint it4 = Back[v3]; it4 < Back[v3 + 1]; ++it4, ++e4) {
          const uint &v4 = e4->idx;
          if (v4 <= st) break;
          if (v4 == v1 || v4 == v2) continue;
          const ulong &w4 = e4->w;
          if (!judge(w4, w3)) continue;
          Data.Reach[v4] |= 8;
          Data.ReachPoint[Data.ReachCount++] = v4;
        }
      }
    }
  }
}

void ForwardSearch1(ThData1 &Data, const uint &st) {
  std::vector<uint>(&ret)[6] = Cycles[st].cycle;

  const DFSEdge *e1 = &G[Head[st]];
  for (uint it1 = Head[st]; it1 < Head[st + 1]; ++it1, ++e1) {
    const uint &v1 = e1->idx;
    const ulong &w1 = e1->w;
    if (v1 < st) continue;

    const DFSEdge *e2 = &G[Head[v1]];

    for (uint it2 = Head[v1]; it2 < Head[v1 + 1]; ++it2, ++e2) {
      const uint &v2 = e2->idx;
      const ulong &w2 = e2->w;
      if (v2 <= st || !judge(w1, w2)) continue;

      const DFSEdge *e3 = &G[Head[v2]];

      for (uint it3 = Head[v2]; it3 < Head[v2 + 1]; ++it3, ++e3) {
        const uint &v3 = e3->idx;
        const ulong &w3 = e3->w;

        if (v3 < st || v3 == v1 || !judge(w2, w3)) {
          continue;
        } else if (v3 == st) {
          if (judge(w3, w1)) {
            ret[0].insert(ret[0].end(), {v1, v2});
          }
          continue;
        }

        const DFSEdge *e4 = &G[Head[v3]];

        for (uint it4 = Head[v3]; it4 < Head[v3 + 1]; ++it4, ++e4) {
          const uint &v4 = e4->idx;
          const ulong &w4 = e4->w;

          if (!(Data.Reach[v4] & 8) || v1 == v4 || v2 == v4) {
            continue;
          } else if (!judge(w3, w4)) {
            continue;
          } else if (v4 == st) {
            if (judge(w4, w1)) {
              ret[1].insert(ret[1].end(), {v1, v2, v3});
            }
            continue;
          }

          const DFSEdge *e5 = &G[Head[v4]];

          for (uint it5 = Head[v4]; it5 < Head[v4 + 1]; ++it5, ++e5) {
            const uint &v5 = e5->idx;
            const ulong &w5 = e5->w;

            if (!(Data.Reach[v5] & 4) || v1 == v5 || v2 == v5 || v3 == v5) {
              continue;
            } else if (!judge(w4, w5)) {
              continue;
            } else if (v5 == st) {
              if (judge(w5, w1)) {
                ret[2].insert(ret[2].end(), {v1, v2, v3, v4});
              }
              continue;
            }

            const DFSEdge *e6 = &G[Head[v5]];

            for (uint it6 = Head[v5]; it6 < Head[v5 + 1]; ++it6, ++e6) {
              const uint &v6 = e6->idx;
              const ulong &w6 = e6->w;

              if (!(Data.Reach[v6] & 2) || v1 == v6 || v2 == v6 || v3 == v6 ||
                  v4 == v6) {
                continue;
              } else if (!judge(w5, w6)) {
                continue;
              } else if (v6 == st) {
                if (judge(w6, w1)) {
                  ret[3].insert(ret[3].end(), {v1, v2, v3, v4, v5});
                }
                continue;
              }

              const DFSEdge *e7 = &G[Head[v6]];

              for (uint it7 = Head[v6]; it7 < Head[v6 + 1]; ++it7, ++e7) {
                const uint &v7 = e7->idx;
                const ulong &w7 = e7->w;

                if (!(Data.Reach[v7] & 1) || v1 == v7 || v2 == v7 || v3 == v7 ||
                    v4 == v7 || v5 == v7) {
                  continue;
                } else if (!judge(w6, w7)) {
                  continue;
                } else if (v7 == st) {
                  if (judge(w7, w1)) {
                    ret[4].insert(ret[4].end(), {v1, v2, v3, v4, v5, v6});
                  }
                  continue;
                }

                const ulong &w8 = Data.LastWeight[v7];

                if (!judge(w7, w8) || !judge(w8, w1)) {
                  continue;
                }

                ret[5].insert(ret[5].end(), {v1, v2, v3, v4, v5, v6, v7});
              }
            }
          }
        }
      }
    }
  }
}

// 8 4 2 1
void BackSearch2(ThData2 &Data, const uint &st) {
  for (uint i = 0; i < Data.ReachCount; ++i) {
    const uint &v = Data.ReachPoint[i];
    Data.Reach[v] = 0;
    Data.Path2[v].resize(0);
    Data.Path3[v].resize(0);
    Data.Path4[v].resize(0);
    Data.vis[v] = false;
  }
  Data.ReachCount = 0;
  Data.ReachPoint[Data.ReachCount++] = st;
  Data.Reach[st] = 15;

  const DFSEdge *e1 = &GBack[Back[st]];
  for (uint it1 = Back[st]; it1 < Back[st + 1]; ++it1, ++e1) {
    const uint &v1 = e1->idx;
    if (v1 <= st) break;
    const ulong &w1 = e1->w;
    Data.LastWeight[v1] = w1;
    Data.Reach[v1] = 15;
    Data.ReachPoint[Data.ReachCount++] = v1;

    const DFSEdge *e2 = &GBack[Back[v1]];
    for (uint it2 = Back[v1]; it2 < Back[v1 + 1]; ++it2, ++e2) {
      const uint &v2 = e2->idx;
      if (v2 <= st) break;
      const ulong &w2 = e2->w;

      if (!judge(w2, w1)) continue;
      Data.Reach[v2] |= 14;
      Data.ReachPoint[Data.ReachCount++] = v2;
      Data.Path2[v2].emplace_back(P2{v1, w2});

      const DFSEdge *e3 = &GBack[Back[v2]];
      for (uint it3 = Back[v2]; it3 < Back[v2 + 1]; ++it3, ++e3) {
        const uint &v3 = e3->idx;
        if (v3 <= st) break;
        if (v3 == v1) continue;
        const ulong &w3 = e3->w;

        if (!judge(w3, w2)) continue;
        Data.Reach[v3] |= 12;
        Data.ReachPoint[Data.ReachCount++] = v3;
        Data.Path3[v3].emplace_back(P3{v2, v1, w3});

        const DFSEdge *e4 = &GBack[Back[v3]];
        for (uint it4 = Back[v3]; it4 < Back[v3 + 1]; ++it4, ++e4) {
          const uint &v4 = e4->idx;
          if (v4 <= st) break;
          if (v4 == v1 || v4 == v2) continue;
          const ulong &w4 = e4->w;
          if (!judge(w4, w3)) continue;
          Data.Reach[v4] |= 8;
          Data.ReachPoint[Data.ReachCount++] = v4;
          Data.Path4[v4].emplace_back(P4{v3, v2, v1, w4});
        }
      }
    }
  }
}

bool cmp2(const P2 &x, const P2 &y) { return x.p1 < y.p1; }
bool cmp3(const P3 &x, const P3 &y) {
  if (x.p1 == y.p1) return x.p2 < y.p2;
  return x.p1 < y.p1;
}
bool cmp4(const P4 &x, const P4 &y) {
  if (x.p1 == y.p1) {
    if (x.p2 == y.p2) return x.p3 < y.p3;
    return x.p2 < y.p2;
  }
  return x.p1 < y.p1;
}

/*
void ForwardSearch2(ThData2 &Data, const uint &st) {
  std::vector<uint>(&ret)[6] = Cycles[st].cycle;

  const DFSEdge *e1 = &G[Head[st]];
  for (uint it1 = Head[st]; it1 < Head[st + 1]; ++it1, ++e1) {
    const uint &v1 = e1->idx;
    const ulong &w1 = e1->w;
    if (v1 < st) continue;

    const DFSEdge *e2 = &G[Head[v1]];

    for (uint it2 = Head[v1]; it2 < Head[v1 + 1]; ++it2, ++e2) {
      const uint &v2 = e2->idx;
      const ulong &w2 = e2->w;
      if (v2 <= st || !judge(w1, w2)) continue;

      const DFSEdge *e3 = &G[Head[v2]];

      for (uint it3 = Head[v2]; it3 < Head[v2 + 1]; ++it3, ++e3) {
        const uint &v3 = e3->idx;
        const ulong &w3 = e3->w;

        if (v3 < st || v3 == v1 || !judge(w2, w3)) {
          continue;
        } else if (v3 == st) {
          if (judge(w3, w1)) {
            ret[0].insert(ret[0].end(), {v1, v2});
          }
          continue;
        }

        const DFSEdge *e4 = &G[Head[v3]];

        for (uint it4 = Head[v3]; it4 < Head[v3 + 1]; ++it4, ++e4) {
          const uint &v4 = e4->idx;
          const ulong &w4 = e4->w;

          if (!(Data.Reach[v4] & 8) || v1 == v4 || v2 == v4) {
            continue;
          } else if (!judge(w3, w4)) {
            continue;
          } else if (v4 == st) {
            if (judge(w4, w1)) {
              ret[1].insert(ret[1].end(), {v1, v2, v3});
            }
            continue;
          }

          const DFSEdge *e5 = &G[Head[v4]];

          for (uint it5 = Head[v4]; it5 < Head[v4 + 1]; ++it5, ++e5) {
            const uint &v5 = e5->idx;
            const ulong &w5 = e5->w;

            if (!(Data.Reach[v5] & 4) || v1 == v5 || v2 == v5 || v3 == v5) {
              continue;
            } else if (!judge(w4, w5)) {
              continue;
            } else if (v5 == st) {
              if (judge(w5, w1)) {
                ret[2].insert(ret[2].end(), {v1, v2, v3, v4});
              }
              continue;
            }

            const DFSEdge *e6 = &G[Head[v5]];

            for (uint it6 = Head[v5]; it6 < Head[v5 + 1]; ++it6, ++e6) {
              const uint &v6 = e6->idx;
              const ulong &w6 = e6->w;

              if (!(Data.Reach[v6] & 2) || v1 == v6 || v2 == v6 || v3 == v6 ||
                  v4 == v6) {
                continue;
              } else if (!judge(w5, w6)) {
                continue;
              } else if (v6 == st) {
                if (judge(w6, w1)) {
                  ret[3].insert(ret[3].end(), {v1, v2, v3, v4, v5});
                }
                continue;
              }

              const DFSEdge *e7 = &G[Head[v6]];

              for (uint it7 = Head[v6]; it7 < Head[v6 + 1]; ++it7, ++e7) {
                const uint &v7 = e7->idx;
                const ulong &w7 = e7->w;

                if (!(Data.Reach[v7] & 1) || v1 == v7 || v2 == v7 || v3 == v7 ||
                    v4 == v7 || v5 == v7) {
                  continue;
                } else if (!judge(w6, w7)) {
                  continue;
                } else if (v7 == st) {
                  if (judge(w7, w1)) {
                    ret[4].insert(ret[4].end(), {v1, v2, v3, v4, v5, v6});
                  }
                  continue;
                }

                const ulong &w8 = Data.LastWeight[v7];

                if (!judge(w7, w8) || !judge(w8, w1)) {
                  continue;
                }

                ret[5].insert(ret[5].end(), {v1, v2, v3, v4, v5, v6, v7});
              }
            }
          }
        }
      }
    }
  }
}
*/
void ForwardSearch2(ThData2 &Data, const uint &st) {
  std::vector<uint>(&ret)[6] = Cycles[st].cycle;

  const DFSEdge *e1 = &G[Head[st]];
  for (uint it1 = Head[st]; it1 < Head[st + 1]; ++it1, ++e1) {
    const uint &v1 = e1->idx;
    const ulong &w1 = e1->w;
    if (v1 < st) continue;

    const DFSEdge *e2 = &G[Head[v1]];

    for (uint it2 = Head[v1]; it2 < Head[v1 + 1]; ++it2, ++e2) {
      const uint &v2 = e2->idx;
      const ulong &w2 = e2->w;
      if (v2 <= st || !judge(w1, w2)) continue;

      const DFSEdge *e3 = &G[Head[v2]];

      for (uint it3 = Head[v2]; it3 < Head[v2 + 1]; ++it3, ++e3) {
        const uint &v3 = e3->idx;
        const ulong &w3 = e3->w;

        if (v3 < st || v3 == v1 || !judge(w2, w3)) {
          continue;
        } else if (v3 == st) {
          if (judge(w3, w1)) {
            ret[0].insert(ret[0].end(), {v1, v2});
          }
          continue;
        }

        const DFSEdge *e4 = &G[Head[v3]];

        for (uint it4 = Head[v3]; it4 < Head[v3 + 1]; ++it4, ++e4) {
          const uint &v4 = e4->idx;
          const ulong &w4 = e4->w;

          if (!(Data.Reach[v4] & 8) || v1 == v4 || v2 == v4) {
            continue;
          } else if (!judge(w3, w4)) {
            continue;
          } else if (v4 == st) {
            if (judge(w4, w1)) {
              ret[1].insert(ret[1].end(), {v1, v2, v3});
            }
            continue;
          }

          const ulong &lastw = Data.LastWeight[v4];
          if ((Data.Reach[v4] & 1) && judge(w4, lastw) && judge(lastw, w1)) {
            ret[2].insert(ret[2].end(), {v1, v2, v3, v4});
          }

          if (!Data.vis[v4]) {
            Data.vis[v4] = true;
            std::sort(Data.Path2[v4].begin(), Data.Path2[v4].end(), cmp2);
            std::sort(Data.Path3[v4].begin(), Data.Path3[v4].end(), cmp3);
            std::sort(Data.Path4[v4].begin(), Data.Path4[v4].end(), cmp4);
          }

          for (const auto &it : Data.Path2[v4]) {
            const uint &v5 = it.p1;
            if (v5 == v1 || v5 == v2 || v5 == v3) continue;
            const ulong &last = Data.LastWeight[v5];
            const ulong &w5 = it.w;
            if (judge(w4, w5) && judge(last, w1)) {
              ret[3].insert(ret[3].end(), {v1, v2, v3, v4, v5});
            }
          }

          for (const auto &it : Data.Path3[v4]) {
            const uint &v5 = it.p1, &v6 = it.p2;
            if (v5 == v1 || v5 == v2 || v5 == v3) continue;
            if (v6 == v1 || v6 == v2 || v6 == v3 || v6 == v4) continue;
            const ulong &w6 = it.w;
            const ulong &last = Data.LastWeight[v6];
            if (judge(w4, w6) && judge(last, w1)) {
              ret[4].insert(ret[4].end(), {v1, v2, v3, v4, v5, v6});
            }
          }

          for (const auto &it : Data.Path4[v4]) {
            const uint &v5 = it.p1, &v6 = it.p2, &v7 = it.p3;
            if (v5 == v1 || v5 == v2 || v5 == v3) continue;
            if (v6 == v1 || v6 == v2 || v6 == v3 || v6 == v4) continue;
            if (v7 == v1 || v7 == v2 || v7 == v3 || v7 == v4 || v7 == v5)
              continue;
            const ulong &w7 = it.w;
            const ulong &last = Data.LastWeight[v7];
            if (judge(w4, w7) && judge(last, w1)) {
              ret[5].insert(ret[5].end(), {v1, v2, v3, v4, v5, v6, v7});
            }
          }
        }
      }
    }
  }
}

// 1: 6+3(剪枝); 2: 4+3(保存路径)
uint ChooseStragety(const uint &job) {
  if (BackLen[job] * 2 > HeadLen[job]) return 1;
  return 2;
  // return random() % 2;
}

void HandleFindCycle(uint pid) {
#ifdef TESTSPEED
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif

  auto &Data1 = ThreadData1[pid];
  auto &Data2 = ThreadData2[pid];
  Data2.Path2 = std::vector<std::vector<P2>>(MaxID);
  Data2.Path3 = std::vector<std::vector<P3>>(MaxID);
  Data2.Path4 = std::vector<std::vector<P4>>(MaxID);

  while (true) {
    while (_JOB_LOCK_.test_and_set())
      ;
    uint job = JobCur < JobsCount ? Jobs[JobCur++] : -1;
    _JOB_LOCK_.clear();
    if (job == -1) break;

    uint choice = ChooseStragety(job);

    if (choice == 1) {
      BackSearch1(Data1, job);
      if (Data1.ReachCount < 2) continue;
      ForwardSearch1(Data1, job);
    } else {
      BackSearch2(Data2, job);
      if (Data2.ReachCount < 2) continue;
      ForwardSearch2(Data2, job);
    }
  }

#ifdef TESTSPEED
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ thread %llu: [find cycle] [cost: %.4fs]\n", pid, t4 - t1);
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
    for (int j = 0; j < 6; ++j) {
      uint sz = ret.cycle[j].size();
      Answers += sz / (j + 2);
      tolSize += sz;
    }
  }

  uint block = tolSize / BUFFERBLOCK + 1;
  uint strow = 0, endrow = 0, stcol = 0, edcol = 0, times = 0;
  for (uint i = 0; i < 6; ++i) {
    for (uint j = 0; j < JobsCount; ++j) {
      if (i == 5 && j == JobsCount - 1) {
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
  OffSet.emplace_back(std::array<uint, 4>{strow, 5, stcol, JobsCount});
#ifdef TESTSPEED
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ CalOffSet: [offset: %lu] [cost: %.4fs]\n", OffSet.size(), t4 - t1);
#endif
}

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
        const auto &st = MapID[Jobs[j]];
        for (const auto &v : cycle) {
          if (idx == 0) {
            memcpy(result, st.str, st.len);
            result += st.len;
          }
          const auto &mpid = MapID[v];
          memcpy(result, mpid.str, mpid.len);
          if (++idx == i + 2) {
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
  printf("@ thread %llu: [memcpy] [cost: %.4fs]\n", pid, t4 - t1);
#endif
}

void WriteAnswer(int pid) {
#ifdef TESTSPEED
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif

  char tmp[NUMLENGTH];
  sprintf(tmp, "%llu\n", Answers);
  FILE *fp = fopen(RESULT, "w");
  fwrite(tmp, 1, strlen(tmp), fp);

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
  uint cnt[6] = {0};
  for (uint i = 0; i < JobsCount; ++i) {
    auto &ret = Cycles[Jobs[i]];
    for (uint j = 0; j < 6; ++j) {
      cnt[j] += ret.cycle[j].size() / (j + 2);
    }
  }
  std::cerr << "@ ans: ";
  for (uint i = 0; i < 6; ++i) std::cerr << cnt[i] << " ";
  std::cerr << "\n";
#endif
  return 0;
}