#include <fcntl.h>
#include <sys/mman.h>
#include <sys/time.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
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
#define TRAIN "./data/18908526/test_data.txt"
#define RESULT "./data/18908526/result.txt"
#else
#define TRAIN "/data/test_data.txt"
#define RESULT "/projects/student/result.txt"
#endif
typedef std::pair<uint, uint> Pair;

/*
 * 常量定义
 */
const uint MAXEDGE = 3000000 + 7;  // 最多边数目
const uint MAXN = MAXEDGE << 1;    // 最多点数目
const uint NTHREAD = 4;            // 线程个数
const uint NUMLENGTH = 11;         // ID最大长度

/*
 * 计数
 */
uint MaxID = 0;       // 最大点
uint Answers = 0;     // 环个数
uint EdgesCount = 0;  // 边数目
uint JobsCount = 0;   // 有效点数目

/*
 * 找环
 */
struct ThData {
  uint answers = 0;          // 环数目
  uint bufsize = 0;          // 字节长度
  uint ReachPointCount = 0;  // 反向可达点数目
  char Reach[MAXN];          // 标记反向可达
  uint ReachPoint[MAXN];     // 可达点集合
  uint LastWeight[MAXN];     // 最后一步权重
} ThreadData[NTHREAD];       // 线程找环

/*
 * 结果
 */
uint TotalBufferSize = 0;                                  // 总buffer大小
uint FirstBufLen = 0;                                      // 换个数bufsize
char FirstBuf[NUMLENGTH];                                  // 环个数buf
std::tuple<uint, uint, uint, uint, uint> Offset[NTHREAD];  // 偏移量
uint CycleBufSize[5][MAXN];             // 每个点每种环sz
std::vector<std::vector<uint>> Cycle0;  // 长度为3的环
std::vector<std::vector<uint>> Cycle1;  // 长度为4的环
std::vector<std::vector<uint>> Cycle2;  // 长度为5的环
std::vector<std::vector<uint>> Cycle3;  // 长度为6的环
std::vector<std::vector<uint>> Cycle4;  // 长度为7的环

/*
 * atomic 锁
 */
uint JobCur = 0;
uint EdgeCur = 0;
uint IDCur = 0;
std::atomic_flag _JOB_LOCK_ = ATOMIC_FLAG_INIT;
std::atomic_flag _EDGE_LOCK_ = ATOMIC_FLAG_INIT;
std::atomic_flag _ID_LOCK_ = ATOMIC_FLAG_INIT;

/*
 * 图信息
 */
struct PreBuffer {
  char str[NUMLENGTH];
  uint len;
};
struct Edge {
  uint u, v, w;
  bool operator<(const Edge &r) const {
    if (u == r.u) return v < r.v;
    return u < r.u;
  }
};
PreBuffer MapID[MAXN];                              // 解析int
uint Jobs[MAXN];                                    // 有效点
uint IDDom[MAXN];                                   // ID集合
Edge Edges[MAXEDGE];                                // 所有边
Pair Children[MAXN];                                // sons
std::vector<Pair> Parents[MAXN];                    // fathers
std::vector<Pair> ThChildren[NTHREAD];              // thsons
std::vector<std::vector<Pair>> ThParents[NTHREAD];  // thfather
std::unordered_map<uint, uint> HashID;              // hashID
std::vector<uint> Rank;                             // rankID

void Init() {
  Cycle0.reserve(MaxID);
  Cycle1.reserve(MaxID);
  Cycle2.reserve(MaxID);
  Cycle3.reserve(MaxID);
  Cycle4.reserve(MaxID);
}

inline void ParseInteger(const uint &x, uint num) {
  auto &mpid = MapID[x];
  sprintf(mpid.str, "%d,", num);
  mpid.len = strlen(mpid.str);
}

inline void GetEdgeId(uint &cur) {
  while (_EDGE_LOCK_.test_and_set())
    ;
  cur = EdgeCur < EdgesCount ? EdgeCur++ : -1;
  _EDGE_LOCK_.clear();
}
void HandleLoadData(uint pid) {
  auto &children = ThChildren[pid];
  auto &parents = ThParents[pid];
  parents.reserve(MaxID);
  children.reserve(MaxID);
  uint cur = 0;
  while (true) {
    GetEdgeId(cur);
    if (cur == -1) break;
    auto &e = Edges[cur];
    uint p1 = Rank[HashID[e.u]];
    uint p2 = Rank[HashID[e.v]];
    e.u = p1, e.v = p2;
    if (children[p1].second == 0) {
      children[p1].first = cur;
    }
    ++children[p1].second;
    parents[p2].emplace_back(std::make_pair(p1, e.w));
  }
}

inline void GetIDId(uint &cur) {
  while (_ID_LOCK_.test_and_set())
    ;
  cur = IDCur < MaxID ? IDCur++ : -1;
  _ID_LOCK_.clear();
}

void HandleCreateGraph() {
  uint cur = 0;
  while (true) {
    GetIDId(cur);
    if (cur == -1) break;
    uint minx = EdgesCount, len = 0;
    for (uint i = 0; i < NTHREAD; ++i) {
      Parents[cur].insert(Parents[cur].end(), ThParents[i][cur].begin(),
                          ThParents[i][cur].end());
      const auto &cdr = ThChildren[i][cur];
      if (cdr.second > 0) {
        minx = std::min(minx, cdr.first);
        len += cdr.second;
      }
    }
    Children[cur].first = minx;
    Children[cur].second = minx + len;
    std::sort(
        Parents[cur].begin(), Parents[cur].end(),
        [](const Pair &e1, const Pair &e2) { return e1.first > e2.first; });
  }
}

void SortAndRank() {
  std::sort(Edges, Edges + EdgesCount);
  std::vector<uint> vec(MaxID);
  for (uint i = 0; i < MaxID; ++i) vec[i] = i;
  std::sort(vec.begin(), vec.end(),
            [&](const uint &x, const uint &y) { return IDDom[x] < IDDom[y]; });
  Rank.reserve(MaxID);
  for (uint i = 0; i < MaxID; ++i) {
    const uint &x = vec[i];
    ParseInteger(i, IDDom[x]);
    Rank[x] = i;
  }
}

void LoadData() {
  uint fd = open(TRAIN, O_RDONLY);
  uint bufsize = lseek(fd, 0, SEEK_END);
  char *buffer = (char *)mmap(NULL, bufsize, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);
  const char *ptr = buffer;

  uint u = 0, v = 0, w = 0;
  HashID.reserve(MAXN);
  while (ptr - buffer < bufsize) {
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
    Edges[EdgesCount].u = u;
    Edges[EdgesCount].v = v;
    Edges[EdgesCount++].w = w;
    if (HashID.find(u) == HashID.end()) {
      IDDom[MaxID] = u;
      HashID.insert({u, MaxID++});
    }
    if (HashID.find(v) == HashID.end()) {
      IDDom[MaxID] = v;
      HashID.insert({v, MaxID++});
    }
    u = v = w = 0;
  }

  SortAndRank();

  // 多线程存图
  std::thread Th[NTHREAD];
  for (uint i = 1; i < NTHREAD; ++i) Th[i] = std::thread(HandleLoadData, i);
  HandleLoadData(0);
  for (uint i = 1; i < NTHREAD; ++i) Th[i].join();
  for (uint i = 1; i < NTHREAD; ++i) Th[i] = std::thread(HandleCreateGraph);
  HandleCreateGraph();
  for (uint i = 1; i < NTHREAD; ++i) Th[i].join();
  for (uint i = 0; i < MaxID; ++i) {
    if (Children[i].second > 0 && !Parents[i].empty()) {
      Jobs[JobsCount++] = i;
    }
  }
  Init();
#ifdef LOCAL
  std::cerr << "@ u: " << MaxID << ", e: " << EdgesCount << "\n";
#endif
}

inline bool judge(const uint &w1, const uint &w2) {
  return (w2 > W5_MAX || w1 <= P5(w2)) && (w1 > W3_MAX || P3(w1) >= w2);
}

void BackSearch(ThData &Data, const uint &st) {
  for (uint i = 0; i < Data.ReachPointCount; ++i) {
    const uint &v = Data.ReachPoint[i];
    Data.Reach[v] = 0;
  }
  Data.ReachPointCount = 0;
  Data.ReachPoint[Data.ReachPointCount++] = st;
  Data.Reach[st] = 7;
  const auto &parent1 = Parents[st];
  for (const auto &it1 : parent1) {
    const uint &v1 = it1.first;
    if (v1 <= st) break;
    const uint &w1 = it1.second;
    Data.LastWeight[v1] = w1;
    Data.Reach[v1] = 7;
    Data.ReachPoint[Data.ReachPointCount++] = v1;
    const auto &parent2 = Parents[v1];
    for (const auto &it2 : parent2) {
      const uint &v2 = it2.first;
      if (v2 <= st) break;
      const uint &w2 = it2.second;
      if (!judge(w2, w1)) continue;
      Data.Reach[v2] |= 6;
      Data.ReachPoint[Data.ReachPointCount++] = v2;
      const auto &parent3 = Parents[v2];
      for (const auto &it3 : parent3) {
        const uint &v3 = it3.first;
        if (v3 <= st) break;
        const uint &w3 = it3.second;
        if (!judge(w3, w2)) continue;
        Data.Reach[v3] |= 4;
        Data.ReachPoint[Data.ReachPointCount++] = v3;
      }
    }
  }
}

void ForwardSearch(ThData &Data, const uint &st) {
  uint ans = 0, sz0 = 0, sz1 = 0, sz2 = 0, sz3 = 0, sz4 = 0;
  const uint &len0 = MapID[st].len;
  auto &cycle0 = Cycle0[st];
  auto &cycle1 = Cycle1[st];
  auto &cycle2 = Cycle2[st];
  auto &cycle3 = Cycle3[st];
  auto &cycle4 = Cycle4[st];
  const auto &cdr1 = Children[st];
  for (uint it1 = cdr1.first; it1 != cdr1.second; ++it1) {
    const auto &e1 = Edges[it1];
    if (e1.v < st) continue;
    const uint &len1 = MapID[e1.v].len;
    const auto &cdr2 = Children[e1.v];
    for (uint it2 = cdr2.first; it2 != cdr2.second; ++it2) {
      const auto &e2 = Edges[it2];
      if (e2.v <= st || !judge(e1.w, e2.w)) continue;
      const uint len = len0 + len1 + MapID[e2.v].len;
      const auto &cdr3 = Children[e2.v];
      for (uint it3 = cdr3.first; it3 != cdr3.second; ++it3) {
        const auto &e3 = Edges[it3];
        if (e3.v < st || e3.v == e1.v || !judge(e2.w, e3.w)) {
          continue;
        } else if (e3.v == st) {
          if (!judge(e3.w, e1.w)) continue;
          cycle0.insert(cycle0.end(), {st, e1.v, e2.v});
          sz0 += len;
          ++ans;
          continue;
        }
        const uint &len3 = MapID[e3.v].len;
        const auto &cdr4 = Children[e3.v];
        for (uint it4 = cdr4.first; it4 != cdr4.second; ++it4) {
          const auto &e4 = Edges[it4];
          if (!(Data.Reach[e4.v] & 4) || !judge(e3.w, e4.w)) {
            continue;
          } else if (e4.v == st) {
            if (!judge(e4.w, e1.w)) continue;
            cycle1.insert(cycle1.end(), {st, e1.v, e2.v, e3.v});
            sz1 += len + len3;
            ++ans;
            continue;
          } else if (e1.v == e4.v || e2.v == e4.v) {
            continue;
          }
          const uint &len4 = MapID[e4.v].len;
          const auto &cdr5 = Children[e4.v];
          for (uint it5 = cdr5.first; it5 != cdr5.second; ++it5) {
            const auto &e5 = Edges[it5];
            if (!(Data.Reach[e5.v] & 2) || !judge(e4.w, e5.w)) {
              continue;
            } else if (e5.v == st) {
              if (!judge(e5.w, e1.w)) continue;
              cycle2.insert(cycle2.end(), {st, e1.v, e2.v, e3.v, e4.v});
              sz2 += len + len3 + len4;
              ++ans;
              continue;
            } else if (e1.v == e5.v || e2.v == e5.v || e3.v == e5.v) {
              continue;
            }
            const uint &len5 = MapID[e5.v].len;
            const auto &cdr6 = Children[e5.v];
            for (uint it6 = cdr6.first; it6 != cdr6.second; ++it6) {
              const auto &e6 = Edges[it6];
              if (!(Data.Reach[e6.v] & 1) || !judge(e5.w, e6.w)) {
                continue;
              } else if (e6.v == st) {
                if (!judge(e6.w, e1.w)) continue;
                cycle3.insert(cycle3.end(), {st, e1.v, e2.v, e3.v, e4.v, e5.v});
                sz3 += len + len3 + len4 + len5;
                ++ans;
                continue;
              }
              const uint &w7 = Data.LastWeight[e6.v];
              if (e1.v == e6.v || e2.v == e6.v || e3.v == e6.v ||
                  e4.v == e6.v || !judge(e6.w, w7) || !judge(w7, e1.w)) {
                continue;
              }
              const uint &len6 = MapID[e6.v].len;
              cycle4.insert(cycle4.end(),
                            {st, e1.v, e2.v, e3.v, e4.v, e5.v, e6.v});
              sz4 += len + len3 + len4 + len5 + len6;
              ++ans;
            }
          }
        }
      }
    }
  }
  Data.answers += ans;
  Data.bufsize += sz0 + sz1 + sz2 + sz3 + sz4;
  CycleBufSize[0][st] = sz0;
  CycleBufSize[1][st] = sz1;
  CycleBufSize[2][st] = sz2;
  CycleBufSize[3][st] = sz3;
  CycleBufSize[4][st] = sz4;
}

inline void GetNextJob(uint &job) {
  while (_JOB_LOCK_.test_and_set())
    ;
  job = JobCur < JobsCount ? Jobs[JobCur++] : -1;
  _JOB_LOCK_.clear();
}

void FindCircle(uint pid) {
#ifdef LOCAL
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

#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ thread %d: Search Cycle, %.4fs\n", pid, t4 - t1);
#endif
}

void CalOffset() {
  uint block = TotalBufferSize / NTHREAD;
  uint tol = 0, x = 0, stl = 0, stidx = 0;
  uint tidx = 0;
  for (uint i = 0; i < 5; ++i) {
    for (uint j = 0; j < JobsCount; ++j) {
      if (x > block) {
        Offset[tidx++] = std::make_tuple(stl, i, stidx, j, tol);
        stl = i;
        stidx = j;
        tol += x;
        x = 0;
      }
      x += CycleBufSize[i][Jobs[j]];
    }
  }
  Offset[tidx++] = std::make_tuple(stl, 4, stidx, JobsCount, tol);
}

void HandleSaveAnswer(uint pid, char *result) {
#ifdef LOCAL
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif

  const auto &offset = Offset[pid];
  uint stl = std::get<0>(offset);
  uint edl = std::get<1>(offset);
  uint stidx = std::get<2>(offset);
  uint edidx = std::get<3>(offset);
  uint offsz = FirstBufLen + std::get<4>(Offset[pid]);
  result += offsz;
  uint low = 0, high = 0;
  auto foo = [&](const uint &p, const std::vector<std::vector<uint>> &cycles) {
    for (uint i = low; i < high; ++i) {
      const auto &cycle = cycles[Jobs[i]];
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
    }
  };
  for (uint i = stl; i <= edl; ++i) {
    low = (i == stl ? stidx : 0);
    high = (i == edl ? edidx : JobsCount);
    switch (i) {
      case 0:
        foo(i, Cycle0);
        break;
      case 1:
        foo(i, Cycle1);
        break;
      case 2:
        foo(i, Cycle2);
        break;
      case 3:
        foo(i, Cycle3);
        break;
      case 4:
        foo(i, Cycle4);
        break;
      default:
        break;
    }
  }

#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ thread %d: Save Answer, %.4fs\n", pid, t4 - t1);
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

void Simulation() {
  LoadData();
  std::thread Th[NTHREAD];
  for (uint i = 1; i < NTHREAD; ++i) Th[i] = std::thread(FindCircle, i);
  FindCircle(0);
  for (uint i = 1; i < NTHREAD; ++i) Th[i].join();

  for (auto &it : ThreadData) {
    Answers += it.answers;
    TotalBufferSize += it.bufsize;
  }

  SaveAnswer();

#ifdef LOCAL
  std::cerr << "@ Answers: " << Answers << "";
  std::cerr << ", Bufsize: " << TotalBufferSize << "\n";
#endif
}

int main() {
  Simulation();
  return 0;
}