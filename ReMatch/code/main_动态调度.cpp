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

#ifdef LOCAL
#define TRAIN "./data/big/test_data.txt"
#define RESULT "./data/big/result.txt"
#else
#define TRAIN "/data/test_data.txt"
#define RESULT "/projects/student/result.txt"
#endif

const int MAXEDGE = 3000000 + 7;
const int MAXNODE = MAXEDGE << 1;
const int NTHREAD = 4;
struct NumBuffer {
  char str[10];
  U32 len;
} MapID[MAXNODE];
struct ThData {
  U32 answers;
  U32 ReachablePointCount;
  char Reachable[MAXNODE];
  U32 ReachablePoint[MAXNODE];
  U32 LastWeight[MAXNODE];
} ThreadData[NTHREAD];

std::atomic_flag lock = ATOMIC_FLAG_INIT;
U32 MaxID = 0;
U32 answers = 0;
U32 jobcur = 0;
U32 IDDomCount = 0;
U32 IDDom[MAXNODE];
std::vector<U32> Jobs;
std::tuple<U32, U32, U32> Edges[MAXEDGE];
U32 EdgesCount;
std::tuple<U32, U32, U32> ThEdges[NTHREAD][MAXEDGE];
U32 ThEdgesCount[NTHREAD];
std::vector<std::pair<U32, U32>> Children[MAXNODE];
std::vector<std::pair<U32, U32>> Parents[MAXNODE];
char *Answer3[MAXNODE];
char *Answer4[MAXNODE];
char *Answer5[MAXNODE];
char *Answer6[MAXNODE];
char *Answer7[MAXNODE];

void ParseInteger(const U32 &x) {
  U32 num = IDDom[x];
  auto &mpid = MapID[x];
  if (num == 0) {
    mpid.str[0] = '0';
    mpid.str[1] = ',';
    mpid.len = 2;
  } else {
    char tmp[10];
    int idx = 10;
    tmp[--idx] = ',';
    while (num) {
      tmp[--idx] = num % 10 + '0';
      num /= 10;
    }
    memcpy(mpid.str, tmp + idx, 10 - idx);
    mpid.len = 10 - idx;
  }
}
void HandleLoadData(int pid, int st, int ed) {
  auto &E = ThEdges[pid];
  auto &count = ThEdgesCount[pid];
  for (int i = st; i < ed; ++i) {
    const auto e = Edges[i];
    const U32 &eu = std::get<0>(e);
    const U32 &ev = std::get<1>(e);
    const U32 &ew = std::get<2>(e);
    auto p1 = std::lower_bound(IDDom, IDDom + MaxID, eu) - IDDom;
    auto p2 = std::lower_bound(IDDom, IDDom + MaxID, ev) - IDDom;
    E[count++] = std::make_tuple(p1, p2, ew);
  }
}
void LoadData() {
  U32 fd = open(TRAIN, O_RDONLY);
  U32 bufsize = lseek(fd, 0, SEEK_END);
  char *buffer = (char *)mmap(NULL, bufsize, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);
  const char *ptr = buffer;
  U32 u = 0, v = 0, w = 0;
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
    while (*ptr != '\r' && *ptr != '\n') {
      w = w * 10 + *ptr - '0';
      ++ptr;
    }
    if (*ptr == '\r') ++ptr;
    ++ptr;
    Edges[EdgesCount++] = std::make_tuple(u, v, w);
    IDDom[IDDomCount++] = u;
    IDDom[IDDomCount++] = v;
    u = v = w = 0;
  }
  std::sort(IDDom, IDDom + IDDomCount);
  MaxID = std::unique(IDDom, IDDom + IDDomCount) - IDDom;
  int st = 0, block = EdgesCount / NTHREAD;
  std::thread Th[NTHREAD];
  for (int i = 0; i < NTHREAD; ++i) {
    int ed = (i == NTHREAD - 1 ? EdgesCount : st + block);
    Th[i] = std::thread(HandleLoadData, i, st, ed);
    st = ed;
  }
  for (int i = 0; i < NTHREAD; ++i) Th[i].join();
  for (int i = 0; i < NTHREAD; ++i) {
    const auto &count = ThEdgesCount[i];
    for (int j = 0; j < count; ++j) {
      const auto &e = ThEdges[i][j];
      const U32 &eu = std::get<0>(e);
      const U32 &ev = std::get<1>(e);
      const U32 &ew = std::get<2>(e);
      Children[eu].emplace_back(std::make_pair(ev, ew));
      Parents[ev].emplace_back(std::make_pair(eu, ew));
    }
  }
  for (int i = 0; i < MaxID; ++i) {
    if (!Children[i].empty() && !Parents[i].empty()) {
      Jobs.emplace_back(i);
      ParseInteger(i);
      std::sort(Children[i].begin(), Children[i].end());
    }
  }
#ifdef LOCAL
  std::cerr << "@ u: " << MaxID << ", e: " << EdgesCount << "\n";
#endif
}

void BackSearch(ThData &Data, const U32 &job) {
  for (int i = 0; i < Data.ReachablePointCount; ++i) {
    Data.Reachable[Data.ReachablePoint[i]] = 0;
  }
  Data.ReachablePointCount = 0;
  Data.ReachablePoint[Data.ReachablePointCount++] = job;
  Data.Reachable[job] = 7;
  U32 v1, v2, v3;
  U32 w1, w2, w3;
  for (const auto &it1 : Parents[job]) {
    v1 = it1.first, w1 = it1.second;
    if (v1 <= job) continue;
    Data.LastWeight[v1] = w1;
    Data.Reachable[v1] = 7;
    Data.ReachablePoint[Data.ReachablePointCount++] = v1;
    for (const auto &it2 : Parents[v1]) {
      v2 = it2.first, w2 = it2.second;
      if (v2 <= job) continue;
      if (w2 > 5 * w1 || 3 * w2 < w1) continue;
      Data.Reachable[v2] |= 6;
      Data.ReachablePoint[Data.ReachablePointCount++] = v2;
      for (const auto &it3 : Parents[v2]) {
        v3 = it3.first, w3 = it3.second;
        if (v3 <= job || v3 == v1) continue;
        if (w3 > 5 * w2 || 3 * w3 < w2) continue;
        Data.Reachable[v3] |= 4;
        Data.ReachablePoint[Data.ReachablePointCount++] = v3;
      }
    }
  }
}
// w2 <= 3w1 && 5w2 >= w1
void ForwardSearch(ThData &Data, const U32 &job) {
  U32 v1, v2, v3, v4, v5, v6;
  U32 w1, w2, w3, w4, w5, w6, w7;
  for (const auto &it1 : Children[job]) {
    v1 = it1.first, w1 = it1.second;
    if (v1 < job) continue;
    for (const auto &it2 : Children[v1]) {
      v2 = it2.first, w2 = it2.second;
      if (v2 <= job || w2 > 3 * w1 || w1 > 5 * w2) continue;
      for (const auto &it3 : Children[v2]) {
        v3 = it3.first, w3 = it3.second;
        if (v3 < job || v3 == v1 || w3 > 3 * w2 || w2 > 5 * w3) continue;
        if (v3 == job) {
          if (w3 <= 5 * w1 && 3 * w3 >= w1) {
            ++Data.answers;
          }
          continue;
        }
        for (const auto &it4 : Children[v3]) {
          v4 = it4.first, w4 = it4.second;
          if (!(Data.Reachable[v4] & 4) || w4 > 3 * w3 || w3 > 5 * w4) continue;
          if (v4 == job) {
            if (w4 <= 5 * w1 && 3 * w4 >= w1) {
              ++Data.answers;
            }
            continue;
          }
          if (v4 == v1 || v4 == v2) {
            continue;
          }
          for (const auto &it5 : Children[v4]) {
            v5 = it5.first, w5 = it5.second;
            if (!(Data.Reachable[v5] & 2) || w5 > 3 * w4 || w4 > 5 * w5)
              continue;
            if (v5 == job) {
              if (w5 <= 5 * w1 && 3 * w5 >= w1) {
                ++Data.answers;
              }
              continue;
            }
            if (v5 == v1 || v5 == v2 || v5 == v3) {
              continue;
            }
            for (const auto &it6 : Children[v5]) {
              v6 = it6.first, w6 = it6.second;
              if (!(Data.Reachable[v6] & 1) || w6 > 3 * w5 || w5 > 5 * w6)
                continue;
              if (v6 == job) {
                if (w6 <= 5 * w1 && 3 * w6 >= w1) {
                  ++Data.answers;
                }
                continue;
              }
              w7 = Data.LastWeight[v6];
              if (v6 == v1 || v6 == v2 || v6 == v3 || v6 == v4 || w7 > 3 * w6 ||
                  w6 > 5 * w7) {
                continue;
              }
              if (w7 <= 5 * w1 && 3 * w7 >= w1) {
                ++Data.answers;
              }
            }
          }
        }
      }
    }
  }
}

void GetNextJob(U32 &job) {
  while (lock.test_and_set())
    ;
  if (jobcur < Jobs.size()) {
    job = Jobs[jobcur++];
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
void Simulation() {
  LoadData();
  std::thread Th[NTHREAD];
  for (int i = 0; i < NTHREAD; ++i) {
    Th[i] = std::thread(FindCircle, i);
  }
  for (auto &it : Th) it.join();
  for (auto &it : ThreadData) answers += it.answers;

#ifdef LOCAL
  std::cerr << "@ answers: " << answers << "\n";
  for (int i = 0; i < NTHREAD; ++i) {
    std::cerr << "@ " << i << ": " << ThreadData[i].answers << "\n";
  }
#endif
}

int main() {
  Simulation();
  return 0;
}