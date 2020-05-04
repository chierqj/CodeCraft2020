#include <fcntl.h>
#include <sys/mman.h>
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
typedef std::pair<uint, uint> Edge;

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

  uint t1 = 0, t2 = 0, t3 = 0;
} ThreadData[NTHREAD];  // 线程找环

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
uint JobCur = 0;  // job光标
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
PreBuffer MapID[MAXN];                               // 解析int
uint Jobs[MAXN];                                     // 有效点
uint IDDom[MAXN];                                    // ID集合
uint Edges[MAXEDGE][3];                              // 所有边
std::vector<Edge> Children[MAXN];                    // 子结点
std::vector<Edge> Parents[MAXN];                     // 父结点
std::vector<std::vector<Edge>> ThChildren[NTHREAD];  // Thread子结点
std::vector<std::vector<Edge>> ThParents[NTHREAD];   // Thread父结点
std::unordered_map<uint, uint> HashID;               // hashID
std::vector<uint> Rank;                              // rankID
std::vector<uint> MaxWeight;                         // MaxWeight
std::vector<uint> MinWeight;                         // MinWeight

void Init() {
  Cycle0.reserve(MaxID);
  Cycle1.reserve(MaxID);
  Cycle2.reserve(MaxID);
  Cycle3.reserve(MaxID);
  Cycle4.reserve(MaxID);
  MaxWeight = std::vector<uint>(MaxID, 0);
  MinWeight = std::vector<uint>(MaxID, std::numeric_limits<uint>::max());
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
  children.reserve(MaxID);
  parents.reserve(MaxID);
  uint cur = 0;
  while (true) {
    GetEdgeId(cur);
    if (cur == -1) break;
    const auto &e = Edges[cur];
    uint p1 = Rank[HashID[e[0]]];
    uint p2 = Rank[HashID[e[1]]];
    children[p1].emplace_back(std::make_pair(p2, e[2]));
    parents[p2].emplace_back(std::make_pair(p1, e[2]));
    MaxWeight[p1] = std::max(MaxWeight[p1], e[2]);
    MinWeight[p1] = std::min(MinWeight[p1], e[2]);
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
    for (uint i = 0; i < NTHREAD; ++i) {
      Children[cur].insert(Children[cur].end(), ThChildren[i][cur].begin(),
                           ThChildren[i][cur].end());
      Parents[cur].insert(Parents[cur].end(), ThParents[i][cur].begin(),
                          ThParents[i][cur].end());
    }
    std::sort(
        Children[cur].begin(), Children[cur].end(),
        [](const Edge &e1, const Edge &e2) { return e1.first < e2.first; });
    std::sort(
        Parents[cur].begin(), Parents[cur].end(),
        [](const Edge &e1, const Edge &e2) { return e1.first > e2.first; });
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
    Edges[EdgesCount][0] = u;
    Edges[EdgesCount][1] = v;
    Edges[EdgesCount++][2] = w;
    if (HashID.find(u) == HashID.end()) {
      IDDom[MaxID] = u;
      HashID[u] = MaxID++;
    }
    if (HashID.find(v) == HashID.end()) {
      IDDom[MaxID] = v;
      HashID[v] = MaxID++;
    }
    u = v = w = 0;
  }

  Init();

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

  // 多线程存图
  std::thread Th[NTHREAD];
  for (uint i = 1; i < NTHREAD; ++i) Th[i] = std::thread(HandleLoadData, i);
  HandleLoadData(0);
  for (uint i = 1; i < NTHREAD; ++i) Th[i].join();
  for (uint i = 1; i < NTHREAD; ++i) Th[i] = std::thread(HandleCreateGraph);
  HandleCreateGraph();
  for (uint i = 1; i < NTHREAD; ++i) Th[i].join();

  for (uint i = 0; i < MaxID; ++i) {
    if (!Children[i].empty() && !Parents[i].empty()) {
      Jobs[JobsCount++] = i;
    }
  }
#ifdef LOCAL
  std::cerr << "@ u: " << MaxID << ", e: " << EdgesCount << "\n";
#endif
}

// w1 <= 5w2 && 3w1 >= w2
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

// 5max < w || min > 3w
// w <= 5max && 3w >= min
inline bool judge1(const uint &v, const uint &w) {
  const uint &maxx = MaxWeight[v];
  const uint &minx = MinWeight[v];
  return (maxx > W5_MAX || w <= P5(maxx)) && (w > W3_MAX || P3(w) >= minx);
}

void ForwardSearch(ThData &Data, const uint &st) {
  uint ans = 0, sz0 = 0, sz1 = 0, sz2 = 0, sz3 = 0, sz4 = 0;
  const uint &len0 = MapID[st].len;
  auto &cycle0 = Cycle0[st];
  auto &cycle1 = Cycle1[st];
  auto &cycle2 = Cycle2[st];
  auto &cycle3 = Cycle3[st];
  auto &cycle4 = Cycle4[st];
  const auto &children1 = Children[st];
  for (const auto &it1 : children1) {
    const uint &v1 = it1.first, &w1 = it1.second;
    if (v1 < st) continue;
    if (!judge1(v1, w1)) {
      // ++Data.t3;
      continue;
    }
    const uint &len1 = MapID[v1].len;
    const auto &children2 = Children[v1];
    for (const auto &it2 : children2) {
      const uint &v2 = it2.first, &w2 = it2.second;
      if (v2 <= st || !judge(w1, w2)) continue;
      if (!judge1(v2, w2)) {
        // ++Data.t3;
        continue;
      }
      const uint len = len0 + len1 + MapID[v2].len;
      const auto &children3 = Children[v2];
      for (const auto &it3 : children3) {
        const uint &v3 = it3.first, &w3 = it3.second;
        if (v3 < st || v3 == v1 || !judge(w2, w3)) {
          continue;
        }
        if (!judge1(v3, w3)) {
          // ++Data.t3;
          continue;
        }
        if (v3 == st) {
          if (!judge(w3, w1)) continue;
          cycle0.insert(cycle0.end(), {st, v1, v2});
          sz0 += len;
          ++ans;
          continue;
        }
        const uint &len3 = MapID[v3].len;
        const auto &children4 = Children[v3];
        for (const auto &it4 : children4) {
          const uint &v4 = it4.first, &w4 = it4.second;
          if (!(Data.Reach[v4] & 4) || !judge(w3, w4)) {
            continue;
          }
          if (!judge1(v4, w4)) {
            // ++Data.t3;
            continue;
          }
          if (v4 == st) {
            if (!judge(w4, w1)) continue;
            cycle1.insert(cycle1.end(), {st, v1, v2, v3});
            sz1 += len + len3;
            ++ans;
            continue;
          } else if (v1 == v4 || v2 == v4) {
            continue;
          }
          const uint &len4 = MapID[v4].len;
          const auto &children5 = Children[v4];
          for (const auto &it5 : children5) {
            const uint &v5 = it5.first, &w5 = it5.second;
            if (!(Data.Reach[v5] & 2) || !judge(w4, w5)) {
              continue;
            }
            if (!judge1(v5, w5)) {
              // ++Data.t3;
              continue;
            }
            if (v5 == st) {
              if (!judge(w5, w1)) continue;
              cycle2.insert(cycle2.end(), {st, v1, v2, v3, v4});
              sz2 += len + len3 + len4;
              ++ans;
              continue;
            } else if (v1 == v5 || v2 == v5 || v3 == v5) {
              continue;
            }
            const uint &len5 = MapID[v5].len;
            const auto &children6 = Children[v5];
            for (const auto &it6 : children6) {
              const uint &v6 = it6.first, &w6 = it6.second;
              if (!(Data.Reach[v6] & 1) || !judge(w5, w6)) {
                continue;
              }
              if (!judge1(v6, w6)) {
                // ++Data.t3;
                continue;
              }
              if (v6 == st) {
                if (!judge(w6, w1)) continue;
                cycle3.insert(cycle3.end(), {st, v1, v2, v3, v4, v5});
                sz3 += len + len3 + len4 + len5;
                ++ans;
                continue;
              }
              const uint &w7 = Data.LastWeight[v6];
              if (v1 == v6 || v2 == v6 || v3 == v6 || v4 == v6 ||
                  !judge(w6, w7) || !judge(w7, w1)) {
                continue;
              }
              const uint &len6 = MapID[v6].len;
              cycle4.insert(cycle4.end(), {st, v1, v2, v3, v4, v5, v6});
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
  uint job = 0;
  auto &Data = ThreadData[pid];
  while (true) {
    GetNextJob(job);
    if (job == -1) break;
    BackSearch(Data, job);
    ForwardSearch(Data, job);
  }
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
#ifdef TEST
  for (uint i = 0; i < NTHREAD; ++i) {
    const auto &e = Offset[i];
    std::cerr << "@ offset: (" << std::get<0>(e) << ", " << std::get<1>(e)
              << ") (" << std::get<2>(e) << ", " << std::get<3>(e) << ") "
              << std::get<4>(e) << "\n";
  }
#endif
}

void HandleSaveAnswer(uint pid) {
  const auto &offset = Offset[pid];
  uint stl = std::get<0>(offset);
  uint edl = std::get<1>(offset);
  uint stidx = std::get<2>(offset);
  uint edidx = std::get<3>(offset);
  uint fd = open(RESULT, O_RDWR | O_CREAT, 0666);
  char *result = (char *)mmap(NULL, TotalBufferSize, PROT_READ | PROT_WRITE,
                              MAP_SHARED, fd, 0);
  ftruncate(fd, TotalBufferSize);
  close(fd);
  if (pid == 0) {
    memcpy(result, FirstBuf, FirstBufLen);
    result += FirstBufLen;
  } else {
    uint offsz = FirstBufLen + std::get<4>(Offset[pid]);
    result += offsz;
  }

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
}

void SaveAnswer() {
  CalOffset();
  sprintf(FirstBuf, "%d\n", Answers);
  FirstBufLen = strlen(FirstBuf);
  TotalBufferSize += FirstBufLen;
  std::thread Th[NTHREAD];
  for (uint i = 1; i < NTHREAD; ++i) Th[i] = std::thread(HandleSaveAnswer, i);
  HandleSaveAnswer(0);
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
#ifdef LOCAL
  uint idx = 0;
  for (auto &it : ThreadData) {
    std::cerr << "@ thread " << idx++ << ": " << it.answers << ", "
              << it.bufsize << ", MaxMin: " << it.t3 << "\n";
  }
  std::cerr << "@ answers: " << Answers << "";
  std::cerr << ", bufsize: " << TotalBufferSize << "\n";
#endif
  SaveAnswer();
}

int main() {
  Simulation();
  return 0;
}