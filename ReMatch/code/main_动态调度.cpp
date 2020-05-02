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
#define U32 uint64_t
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

/*
 * 常量定义
 */
const U32 MAXEDGE = 3000000 + 7;  // 最多边数目
const U32 MAXN = MAXEDGE << 1;    // 最多点数目
const int NTHREAD = 4;            // 线程个数
const int NUMLENGTH = 12;         // ID最大长度

/*
 * 计数
 */
U32 MaxID = 0;              // 最大点
U32 Answers = 0;            // 环个数
U32 Jobcur = 0;             // job光标
U32 EdgesCount = 0;         // 边数目
U32 JobsCount = 0;          // 有效点数目
U32 IDDomCount = 0;         // ID数目
U32 ThEdgesCount[NTHREAD];  // 线程边数目

/*
 * 图信息
 */
struct PreBuffer {
  char str[NUMLENGTH];
  U32 len;
};
std::vector<PreBuffer> MapID;                            // 解析int
std::vector<U32> Jobs;                                   // 有效点
U32 IDDom[MAXN];                                         // ID集合
U32 Edges[MAXEDGE][3];                                   // 所有边
U32 ThEdges[NTHREAD][MAXEDGE / NTHREAD + 7][3];          // 线程边
std::vector<std::vector<std::pair<U32, U32>>> Children;  // 子结点
std::vector<std::vector<std::pair<U32, U32>>> Parents;   // 子结点
// U32 Jobs[MAXN];                                          // 有效点

/*
 * 找环
 */
struct ThData {
  U32 answers = 0;                  // 环数目
  U32 bufsize = 0;                  // 字节长度
  U32 ReachablePointCount = 0;      // 反向可达点数目
  char Reachable[MAXN];             // 标记反向可达
  std::vector<U32> ReachablePoint;  // 可达点集合
  std::vector<U32> LastWeight;      // 最后一步权重
} ThreadData[NTHREAD];              // 线程找环

/*
 * 结果
 */
U32 TotalBufferSize = 0;                              // 总buffer大小
U32 FirstBufLen = 0;                                  // 换个数bufsize
char FirstBuf[NUMLENGTH];                             // 环个数buf
std::tuple<int, int, int, int, U32> Offset[NTHREAD];  // 偏移量
std::vector<U32> CycleBufSize[5];                     // 每个点每种环sz
std::vector<std::vector<U32>> Cycle0;                 // 长度为3的环
std::vector<std::vector<U32>> Cycle1;                 // 长度为4的环
std::vector<std::vector<U32>> Cycle2;                 // 长度为5的环
std::vector<std::vector<U32>> Cycle3;                 // 长度为6的环
std::vector<std::vector<U32>> Cycle4;                 // 长度为7的环

/*
 * atomic 锁
 */
std::atomic_flag lock = ATOMIC_FLAG_INIT;

void Init() {
  Children.reserve(MaxID);
  Parents.reserve(MaxID);
  MapID.reserve(MaxID);
  Jobs.reserve(MaxID);
  for (auto &it : ThreadData) {
    it.ReachablePoint.reserve(MaxID);
    it.LastWeight.reserve(MaxID);
  }
  Cycle0.reserve(MaxID);
  Cycle1.reserve(MaxID);
  Cycle2.reserve(MaxID);
  Cycle3.reserve(MaxID);
  Cycle4.reserve(MaxID);
  for (int i = 0; i < 5; ++i) {
    CycleBufSize[i].reserve(MaxID);
  }
}

void ParseInteger(const U32 &x) {
  U32 num = IDDom[x];
  auto &mpid = MapID[x];
  if (num == 0) {
    mpid.str[0] = '0';
    mpid.str[1] = ',';
    mpid.len = 2;
  } else {
    char tmp[NUMLENGTH];
    int idx = NUMLENGTH;
    tmp[--idx] = ',';
    while (num) {
      tmp[--idx] = num % 10 + '0';
      num /= 10;
    }
    memcpy(mpid.str, tmp + idx, NUMLENGTH - idx);
    mpid.len = NUMLENGTH - idx;
  }
}
void HandleLoadData(int pid, int st, int ed) {
  auto &E = ThEdges[pid];
  auto &count = ThEdgesCount[pid];
  for (int i = st; i < ed; ++i) {
    const auto &e = Edges[i];
    E[count][0] = std::lower_bound(IDDom, IDDom + MaxID, e[0]) - IDDom;
    E[count][1] = std::lower_bound(IDDom, IDDom + MaxID, e[1]) - IDDom;
    E[count++][2] = e[2];
  }
}

void LoadData() {
  std::cerr << "loadstart\n";
  U32 fd = open(TRAIN, O_RDONLY);
  U32 bufsize = lseek(fd, 0, SEEK_END);
  char *buffer = (char *)mmap(NULL, bufsize, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);
  const char *ptr = buffer;
  U32 u = 0, v = 0, w = 0;
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
    IDDom[IDDomCount++] = u;
    IDDom[IDDomCount++] = v;
    u = v = w = 0;
  }
  std::sort(IDDom, IDDom + IDDomCount);
  MaxID = std::unique(IDDom, IDDom + IDDomCount) - IDDom;
  Init();
  int st = 0, block = EdgesCount / NTHREAD;
  std::thread Th[NTHREAD];
  for (int i = 0; i < NTHREAD; ++i) {
    int ed = (i == NTHREAD - 1 ? EdgesCount : st + block);
    Th[i] = std::thread(HandleLoadData, i, st, ed);
    st = ed;
  }
  for (int i = 0; i < NTHREAD; ++i) Th[i].join();
  for (int i = 0; i < NTHREAD; ++i) {
    const auto &E = ThEdges[i];
    const auto &count = ThEdgesCount[i];
    for (int j = 0; j < count; ++j) {
      const auto &e = E[j];
      Children[e[0]].emplace_back(std::make_pair(e[1], e[2]));
      Parents[e[1]].emplace_back(std::make_pair(e[0], e[2]));
    }
  }
  for (int i = 0; i < MaxID; ++i) {
    if (!Children[i].empty() && !Parents[i].empty()) {
      Jobs[JobsCount++] = i;
      ParseInteger(i);
      std::sort(Children[i].begin(), Children[i].end());
    }
  }
#ifdef LOCAL
  std::cerr << "@ u: " << MaxID << ", e: " << EdgesCount << "\n";
#endif
}

inline bool judge(const U32 &w1, const U32 &w2) {
  if (w2 > P3(w1) || P5(w2) < w1) return false;
  return true;
}

void BackSearch(ThData &Data, const U32 &st) {
  for (int i = 0; i < Data.ReachablePointCount; ++i) {
    Data.Reachable[Data.ReachablePoint[i]] = 0;
  }
  Data.ReachablePointCount = 0;
  Data.ReachablePoint[Data.ReachablePointCount++] = st;
  Data.Reachable[st] = 7;
  for (const auto &it1 : Parents[st]) {
    const U32 &v1 = it1.first, &w1 = it1.second;
    if (v1 <= st) continue;
    Data.LastWeight[v1] = w1;
    Data.Reachable[v1] = 7;
    Data.ReachablePoint[Data.ReachablePointCount++] = v1;
    for (const auto &it2 : Parents[v1]) {
      const U32 &v2 = it2.first, &w2 = it2.second;
      if (v2 <= st || !judge(w2, w1)) continue;
      Data.Reachable[v2] |= 6;
      Data.ReachablePoint[Data.ReachablePointCount++] = v2;
      for (const auto &it3 : Parents[v2]) {
        const U32 &v3 = it3.first, &w3 = it3.second;
        if (v3 <= st || v3 == v1 || !judge(w3, w2)) continue;
        Data.Reachable[v3] |= 4;
        Data.ReachablePoint[Data.ReachablePointCount++] = v3;
      }
    }
  }
}

void ForwardSearch(ThData &Data, const U32 &st) {
  U32 ans = 0, sz0 = 0, sz1 = 0, sz2 = 0, sz3 = 0, sz4 = 0;
  const U32 &len0 = MapID[st].len;
  const auto &mpid0 = MapID[st];
  for (const auto &it1 : Children[st]) {
    const U32 &v1 = it1.first, &w1 = it1.second;
    if (v1 < st) continue;
    const U32 &len1 = MapID[v1].len;
    for (const auto &it2 : Children[v1]) {
      const U32 &v2 = it2.first, &w2 = it2.second;
      if (v2 <= st || !judge(w1, w2)) continue;
      const U32 &len2 = MapID[v2].len;
      for (const auto &it3 : Children[v2]) {
        const U32 &v3 = it3.first, &w3 = it3.second;
        if (v3 < st || v3 == v1 || !judge(w2, w3)) {
          continue;
        } else if (v3 == st) {
          if (!judge(w3, w1)) continue;
          Cycle0[st].insert(Cycle0[st].end(), {st, v1, v2});
          sz0 += len0 + len1 + len2;
          ++ans;
          continue;
        }
        const U32 &len3 = MapID[v3].len;
        for (const auto &it4 : Children[v3]) {
          const U32 &v4 = it4.first, &w4 = it4.second;
          if (!(Data.Reachable[v4] & 4) || !judge(w3, w4)) {
            continue;
          } else if (v4 == st) {
            if (!judge(w4, w1)) continue;
            Cycle1[st].insert(Cycle1[st].end(), {st, v1, v2, v3});
            sz1 += len0 + len1 + len2 + len3;
            ++ans;
            continue;
          } else if (v1 == v4 || v2 == v4) {
            continue;
          }
          const U32 &len4 = MapID[v4].len;
          for (const auto &it5 : Children[v4]) {
            const U32 &v5 = it5.first, &w5 = it5.second;
            if (!(Data.Reachable[v5] & 2) || !judge(w4, w5)) {
              continue;
            } else if (v5 == st) {
              if (!judge(w5, w1)) continue;
              Cycle2[st].insert(Cycle2[st].end(), {st, v1, v2, v3, v4});
              sz2 += len0 + len1 + len2 + len3 + len4;
              ++ans;
              continue;
            } else if (v1 == v5 || v2 == v5 || v3 == v5) {
              continue;
            }
            const U32 &len5 = MapID[v5].len;
            for (const auto &it6 : Children[v5]) {
              const U32 &v6 = it6.first, &w6 = it6.second;
              if (!(Data.Reachable[v6] & 1) || !judge(w5, w6)) {
                continue;
              } else if (v6 == st) {
                if (!judge(w6, w1)) continue;
                Cycle3[st].insert(Cycle3[st].end(), {st, v1, v2, v3, v4, v5});
                sz3 += len0 + len1 + len2 + len3 + len4 + len5;
                ++ans;
                continue;
              }
              const U32 &w7 = Data.LastWeight[v6];
              if (v1 == v6 || v2 == v6 || v3 == v6 || v4 == v6 ||
                  !judge(w6, w7) || !judge(w7, w1)) {
                continue;
              }
              const U32 &len6 = MapID[v6].len;
              Cycle4[st].insert(Cycle4[st].end(), {st, v1, v2, v3, v4, v5, v6});
              sz4 += len0 + len1 + len2 + len3 + len4 + len5 + len6;
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

void GetNextJob(U32 &job) {
  while (lock.test_and_set())
    ;
  if (Jobcur < JobsCount) {
    job = Jobs[Jobcur++];
  } else {
    job = -1;
  }
  lock.clear();
}

void FindCircle(int pid) {
  U32 job = 0;
  auto &Data = ThreadData[pid];
  while (true) {
    GetNextJob(job);
    if (job == -1) break;
    BackSearch(Data, job);
    ForwardSearch(Data, job);
  }
}

void CalOffset() {
  U32 block = TotalBufferSize / NTHREAD;
  U32 tol = 0, x = 0, stl = 0, stidx = 0;
  int tidx = 0;
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < JobsCount; ++j) {
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
#ifdef LOCAL
  for (int i = 0; i < NTHREAD; ++i) {
    const auto &e = Offset[i];
    std::cerr << "@ offset: (" << std::get<0>(e) << ", " << std::get<1>(e)
              << ") (" << std::get<2>(e) << ", " << std::get<3>(e) << ") "
              << std::get<4>(e) << "\n";
  }
#endif
}

void HandleSaveAnswer(int pid) {
  const auto &offset = Offset[pid];
  int stl = std::get<0>(offset);
  int edl = std::get<1>(offset);
  int stidx = std::get<2>(offset);
  int edidx = std::get<3>(offset);
  int fd = open(RESULT, O_RDWR | O_CREAT, 0666);
  char *result = (char *)mmap(NULL, TotalBufferSize, PROT_READ | PROT_WRITE,
                              MAP_SHARED, fd, 0);
  ftruncate(fd, TotalBufferSize);
  close(fd);
  if (pid == 0) {
    memcpy(result, FirstBuf, FirstBufLen);
    result += FirstBufLen;
  } else {
    U32 offsz = FirstBufLen + std::get<4>(Offset[pid]);
    result += offsz;
  }
  for (int i = stl; i <= edl; ++i) {
    int low = (i == stl ? stidx : 0);
    int high = (i == edl ? edidx : JobsCount);
    if (i == 0) {
      for (int j = low; j < high; ++j) {
        const auto &job = Jobs[j];
        int idx = 0;
        for (auto &v : Cycle0[job]) {
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
      for (int j = low; j < high; ++j) {
        const auto &job = Jobs[j];
        int idx = 0;
        for (auto &v : Cycle1[job]) {
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
      for (int j = low; j < high; ++j) {
        const auto &job = Jobs[j];
        int idx = 0;
        for (auto &v : Cycle2[job]) {
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
      for (int j = low; j < high; ++j) {
        const auto &job = Jobs[j];
        int idx = 0;
        for (auto &v : Cycle3[job]) {
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
      for (int j = low; j < high; ++j) {
        const auto &job = Jobs[j];
        int idx = 0;
        for (auto &v : Cycle4[job]) {
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
  U32 idx = NUMLENGTH;
  tmp[--idx] = '\n';
  U32 x = Answers;
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
  for (int i = 0; i < NTHREAD; ++i) {
    Th[i] = std::thread(HandleSaveAnswer, i);
  }
  for (int i = 0; i < NTHREAD; ++i) Th[i].join();
}

void Simulation() {
  LoadData();
  std::thread Th[NTHREAD];
  for (int i = 0; i < NTHREAD; ++i) {
    Th[i] = std::thread(FindCircle, i);
  }
  for (auto &it : Th) it.join();
  for (auto &it : ThreadData) {
    Answers += it.answers;
    TotalBufferSize += it.bufsize;
  }
#ifdef LOCAL
  std::cerr << "@ answers: " << Answers << "";
  std::cerr << ", bufsize: " << TotalBufferSize << "\n";
  int idx = 0;
  for (auto &it : ThreadData) {
    std::cerr << "@ thread " << idx++ << ": " << it.answers << ", "
              << it.bufsize << "\n";
  }
#endif
  SaveAnswer();
}

int main() {
  Simulation();
  return 0;
}