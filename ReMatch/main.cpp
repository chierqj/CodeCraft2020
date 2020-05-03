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
#define ulong uint64_t
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
typedef std::pair<uint, ulong> Edge;

/*
 * 常量定义
 */
const uint MAXEDGE = 3000000 + 7;  // 最多边数目
const uint MAXN = MAXEDGE << 1;    // 最多点数目
const uint NTHREAD = 4;            // 线程个数
const uint NUMLENGTH = 12;         // ID最大长度

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
  ulong LastWeight[MAXN];    // 最后一步权重
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

void Init() {
  Cycle0.reserve(MaxID);
  Cycle1.reserve(MaxID);
  Cycle2.reserve(MaxID);
  Cycle3.reserve(MaxID);
  Cycle4.reserve(MaxID);
  Rank.reserve(MaxID);
}

void ParseInteger(const uint &x, uint num) {
  auto &mpid = MapID[x];
  if (num == 0) {
    mpid.str[0] = '0';
    mpid.str[1] = ',';
    mpid.len = 2;
  } else {
    char tmp[NUMLENGTH];
    uint idx = NUMLENGTH;
    tmp[--idx] = ',';
    while (num) {
      tmp[--idx] = num % 10 + '0';
      num /= 10;
    }
    memcpy(mpid.str, tmp + idx, NUMLENGTH - idx);
    mpid.len = NUMLENGTH - idx;
  }
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
  for (uint i = 0; i < MaxID; ++i) {
    const uint &x = vec[i];
    ParseInteger(i, IDDom[x]);
    Rank[x] = i;
  }

  // 多线程存图
  std::thread Th[NTHREAD];
  for (uint i = 0; i < NTHREAD; ++i) Th[i] = std::thread(HandleLoadData, i);
  for (uint i = 0; i < NTHREAD; ++i) Th[i].join();
  for (uint i = 0; i < NTHREAD; ++i) Th[i] = std::thread(HandleCreateGraph);
  for (uint i = 0; i < NTHREAD; ++i) Th[i].join();

  for (uint i = 0; i < MaxID; ++i) {
    if (!Children[i].empty() && !Parents[i].empty()) {
      Jobs[JobsCount++] = i;
      std::sort(Children[i].begin(), Children[i].end());
    }
  }
#ifdef LOCAL
  std::cerr << "@ u: " << MaxID << ", e: " << EdgesCount << "\n";
#endif
}

inline bool judge(const ulong &w1, const ulong &w2) {
  if (w2 > P3(w1) || P5(w2) < w1) return false;
  return true;
}

void BackSearch(ThData &Data, const uint &st) {
  for (uint i = 0; i < Data.ReachPointCount; ++i) {
    const uint &v = Data.ReachPoint[i];
    Data.Reach[v] = 0;
  }
  Data.ReachPointCount = 0;
  Data.ReachPoint[Data.ReachPointCount++] = st;
  Data.Reach[st] = 7;
  for (const auto &it1 : Parents[st]) {
    const uint &v1 = it1.first;
    const ulong &w1 = it1.second;
    if (v1 <= st) continue;
    Data.LastWeight[v1] = w1;
    Data.Reach[v1] = 7;
    Data.ReachPoint[Data.ReachPointCount++] = v1;
    for (const auto &it2 : Parents[v1]) {
      const uint &v2 = it2.first;
      const ulong &w2 = it2.second;
      if (v2 <= st || !judge(w2, w1)) continue;
      Data.Reach[v2] |= 6;
      Data.ReachPoint[Data.ReachPointCount++] = v2;
      for (const auto &it3 : Parents[v2]) {
        const uint &v3 = it3.first;
        const ulong &w3 = it3.second;
        if (v3 <= st || v3 == v1 || !judge(w3, w2)) continue;
        Data.Reach[v3] |= 4;
        Data.ReachPoint[Data.ReachPointCount++] = v3;
      }
    }
  }
}

void ForwardSearch(ThData &Data, const uint &st) {
  uint ans = 0, sz0 = 0, sz1 = 0, sz2 = 0, sz3 = 0, sz4 = 0;
  const uint &len0 = MapID[st].len;
  for (const auto &it1 : Children[st]) {
    const uint &v1 = it1.first;
    const ulong &w1 = it1.second;
    if (v1 < st) continue;
    const uint &len1 = MapID[v1].len;
    for (const auto &it2 : Children[v1]) {
      const uint &v2 = it2.first;
      const ulong &w2 = it2.second;
      if (v2 <= st || !judge(w1, w2)) continue;
      const uint len = len0 + len1 + MapID[v2].len;
      for (const auto &it3 : Children[v2]) {
        const uint &v3 = it3.first;
        const ulong &w3 = it3.second;
        if (v3 < st || v3 == v1 || !judge(w2, w3)) {
          continue;
        } else if (v3 == st) {
          if (!judge(w3, w1)) continue;
          Cycle0[st].insert(Cycle0[st].end(), {st, v1, v2});
          sz0 += len;
          ++ans;
          continue;
        }
        const uint &len3 = MapID[v3].len;
        for (const auto &it4 : Children[v3]) {
          const uint &v4 = it4.first;
          const ulong &w4 = it4.second;
          if (!(Data.Reach[v4] & 4) || !judge(w3, w4)) {
            continue;
          } else if (v4 == st) {
            if (!judge(w4, w1)) continue;
            Cycle1[st].insert(Cycle1[st].end(), {st, v1, v2, v3});
            sz1 += len + len3;
            ++ans;
            continue;
          } else if (v1 == v4 || v2 == v4) {
            continue;
          }
          const uint &len4 = MapID[v4].len;
          for (const auto &it5 : Children[v4]) {
            const uint &v5 = it5.first;
            const ulong &w5 = it5.second;
            if (!(Data.Reach[v5] & 2) || !judge(w4, w5)) {
              continue;
            } else if (v5 == st) {
              if (!judge(w5, w1)) continue;
              Cycle2[st].insert(Cycle2[st].end(), {st, v1, v2, v3, v4});
              sz2 += len + len3 + len4;
              ++ans;
              continue;
            } else if (v1 == v5 || v2 == v5 || v3 == v5) {
              continue;
            }
            const uint &len5 = MapID[v5].len;
            for (const auto &it6 : Children[v5]) {
              const uint &v6 = it6.first;
              const ulong &w6 = it6.second;
              if (!(Data.Reach[v6] & 1) || !judge(w5, w6)) {
                continue;
              } else if (v6 == st) {
                if (!judge(w6, w1)) continue;
                Cycle3[st].insert(Cycle3[st].end(), {st, v1, v2, v3, v4, v5});
                sz3 += len + len3 + len4 + len5;
                ++ans;
                continue;
              }
              const ulong &w7 = Data.LastWeight[v6];
              if (v1 == v6 || v2 == v6 || v3 == v6 || v4 == v6 ||
                  !judge(w6, w7) || !judge(w7, w1)) {
                continue;
              }
              const uint &len6 = MapID[v6].len;
              Cycle4[st].insert(Cycle4[st].end(), {st, v1, v2, v3, v4, v5, v6});
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
  for (uint i = stl; i <= edl; ++i) {
    uint low = (i == stl ? stidx : 0);
    uint high = (i == edl ? edidx : JobsCount);
    if (i == 0) {
      for (uint j = low; j < high; ++j) {
        const auto &job = Jobs[j];
        uint idx = 0;
        for (const auto &v : Cycle0[job]) {
          ++idx;
          const auto &mpid = MapID[v];
          memcpy(result, mpid.str, mpid.len);
          if (idx == i + 3) {
            idx = 0;
            *(result + mpid.len - 1) = '\n';
          }
          result += mpid.len;
        }
      }
    } else if (i == 1) {
      for (uint j = low; j < high; ++j) {
        const auto &job = Jobs[j];
        uint idx = 0;
        for (const auto &v : Cycle1[job]) {
          ++idx;
          const auto &mpid = MapID[v];
          memcpy(result, mpid.str, mpid.len);
          if (idx == i + 3) {
            idx = 0;
            *(result + mpid.len - 1) = '\n';
          }
          result += mpid.len;
        }
      }
    } else if (i == 2) {
      for (uint j = low; j < high; ++j) {
        const auto &job = Jobs[j];
        uint idx = 0;
        for (const auto &v : Cycle2[job]) {
          ++idx;
          const auto &mpid = MapID[v];
          memcpy(result, mpid.str, mpid.len);
          if (idx == i + 3) {
            idx = 0;
            *(result + mpid.len - 1) = '\n';
          }
          result += mpid.len;
        }
      }
    } else if (i == 3) {
      for (uint j = low; j < high; ++j) {
        const auto &job = Jobs[j];
        uint idx = 0;
        for (const auto &v : Cycle3[job]) {
          ++idx;
          const auto &mpid = MapID[v];
          memcpy(result, mpid.str, mpid.len);
          if (idx == i + 3) {
            idx = 0;
            *(result + mpid.len - 1) = '\n';
          }
          result += mpid.len;
        }
      }
    } else {
      for (uint j = low; j < high; ++j) {
        const auto &job = Jobs[j];
        uint idx = 0;
        for (const auto &v : Cycle4[job]) {
          ++idx;
          const auto &mpid = MapID[v];
          memcpy(result, mpid.str, mpid.len);
          if (idx == i + 3) {
            idx = 0;
            *(result + mpid.len - 1) = '\n';
          }
          result += mpid.len;
        }
      }
    }
  }
}

void SaveAnswer() {
  CalOffset();
  char tmp[NUMLENGTH];
  uint idx = NUMLENGTH;
  tmp[--idx] = '\n';
  uint x = Answers;
  if (x == 0) {
    tmp[--idx] = '0';
  } else {
    while (x) {
      tmp[--idx] = x % 10 + '0';
      x /= 10;
    }
  }
  FirstBufLen = NUMLENGTH - idx;
  memcpy(FirstBuf, tmp + idx, FirstBufLen);
  TotalBufferSize += FirstBufLen;
  std::thread Th[NTHREAD];
  for (uint i = 0; i < NTHREAD; ++i) {
    Th[i] = std::thread(HandleSaveAnswer, i);
  }
  for (uint i = 0; i < NTHREAD; ++i) Th[i].join();
}

void Simulation() {
  LoadData();

  std::thread Th[NTHREAD];
  for (uint i = 0; i < NTHREAD; ++i) Th[i] = std::thread(FindCircle, i);
  for (auto &it : Th) it.join();

  for (auto &it : ThreadData) {
    Answers += it.answers;
    TotalBufferSize += it.bufsize;
  }
#ifdef LOCAL
  uint idx = 0;
  for (auto &it : ThreadData) {
    std::cerr << "@ thread " << idx++ << ": " << it.answers << ", "
              << it.bufsize << "\n";
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