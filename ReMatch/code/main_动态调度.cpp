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
const U32 MAXEDGE = 3000000 + 7;    // 最多边数目
const U32 MAXN = MAXEDGE << 1;      // 最多点数目
const U32 MAXCYCLE = 20000000 + 7;  // 最多环数
const int NTHREAD = 4;              // 线程个数
const int NUMLENGTH = 12;           // ID最大长度

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
} MapID[MAXN];                                    // 解析int
U32 Jobs[MAXN];                                   // 有效点
U32 IDDom[MAXN];                                  // ID集合
U32 Edges[MAXEDGE][3];                            // 所有边
U32 ThEdges[NTHREAD][MAXEDGE / NTHREAD + 7][3];   // 线程边
std::vector<std::pair<U32, U32>> Children[MAXN];  // 子结点
std::vector<std::pair<U32, U32>> Parents[MAXN];   // 父节点

/*
 * 找环
 */
struct ThData {
  U32 answers;                                        // 环数目
  U32 ReachablePointCount;                            // 反向可达点数目
  char Reachable[MAXN];                               // 标记反向可达
  U32 ReachablePoint[MAXN];                           // 可达点集合
  U32 LastWeight[MAXN];                               // 最后一步权重
  char *Ans0, *Ans1, *Ans2, *Ans3, *Ans4;             // ans
  char *ThreeCycle = new char[NUMLENGTH * 3];         // 长度为3的环
  char *Cycle0 = new char[MAXCYCLE * NUMLENGTH * 3];  // 长度为3的环
  char *Cycle1 = new char[MAXCYCLE * NUMLENGTH * 4];  // 长度为4的环
  char *Cycle2 = new char[MAXCYCLE * NUMLENGTH * 5];  // 长度为5的环
  char *Cycle3 = new char[MAXCYCLE * NUMLENGTH * 6];  // 长度为6的环
  char *Cycle4 = new char[MAXCYCLE * NUMLENGTH * 7];  // 长度为7的环
} ThreadData[NTHREAD];                                // 线程找环

/*
 * atomic 锁
 */
std::atomic_flag lock = ATOMIC_FLAG_INIT;

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
  U32 ans = 0;
  for (const auto &it1 : Children[st]) {
    const U32 &v1 = it1.first, &w1 = it1.second;
    if (v1 < st) continue;
    for (const auto &it2 : Children[v1]) {
      const U32 &v2 = it2.first, &w2 = it2.second;
      if (v2 <= st || !judge(w1, w2)) continue;
      const auto &mpid0 = MapID[st];
      const auto &mpid1 = MapID[v1];
      const auto &mpid2 = MapID[v2];
      memcpy(Data.ThreeCycle, mpid0.str, mpid0.len);
      memcpy(Data.ThreeCycle + mpid0.len, mpid1.str, mpid1.len);
      memcpy(Data.ThreeCycle + mpid0.len + mpid1.len, mpid2.str, mpid2.len);
      U32 ThreeLength = mpid0.len + mpid1.len + mpid2.len;
      for (const auto &it3 : Children[v2]) {
        const U32 &v3 = it3.first, &w3 = it3.second;
        if (v3 < st || v3 == v1 || !judge(w2, w3)) {
          continue;
        } else if (v3 == st) {
          if (!judge(w3, w1)) continue;
          memcpy(Data.Ans0, Data.ThreeCycle, ThreeLength);
          *(Data.Ans0 + ThreeLength - 1) = '\n';
          Data.Ans0 += ThreeLength;
          ++ans;
          continue;
        }
        const auto &mpid3 = MapID[v3];
        for (const auto &it4 : Children[v3]) {
          const U32 &v4 = it4.first, &w4 = it4.second;
          if (!(Data.Reachable[v4] & 4) || !judge(w3, w4)) {
            continue;
          } else if (v4 == st) {
            if (!judge(w4, w1)) continue;
            memcpy(Data.Ans1, Data.ThreeCycle, ThreeLength);
            Data.Ans1 += ThreeLength;
            memcpy(Data.Ans1, mpid3.str, mpid3.len);
            *(Data.Ans1 + mpid3.len - 1) = '\n';
            Data.Ans1 += mpid3.len;
            ++ans;
            continue;
          } else if (v1 == v4 || v2 == v4) {
            continue;
          }
          const auto &mpid4 = MapID[v4];
          for (const auto &it5 : Children[v4]) {
            const U32 &v5 = it5.first, &w5 = it5.second;
            if (!(Data.Reachable[v5] & 2) || !judge(w4, w5)) {
              continue;
            } else if (v5 == st) {
              if (!judge(w5, w1)) continue;
              memcpy(Data.Ans2, Data.ThreeCycle, ThreeLength);
              Data.Ans2 += ThreeLength;
              memcpy(Data.Ans2, mpid3.str, mpid3.len);
              Data.Ans2 += mpid3.len;
              memcpy(Data.Ans2, mpid4.str, mpid4.len);
              *(Data.Ans2 + mpid4.len - 1) = '\n';
              Data.Ans2 += mpid4.len;
              ++ans;
              continue;
            } else if (v1 == v5 || v2 == v5 || v3 == v5) {
              continue;
            }
            const auto &mpid5 = MapID[v5];
            for (const auto &it6 : Children[v5]) {
              const U32 &v6 = it6.first, &w6 = it6.second;
              if (!(Data.Reachable[v6] & 1) || !judge(w5, w6)) {
                continue;
              } else if (v6 == st) {
                if (!judge(w6, w1)) continue;
                memcpy(Data.Ans3, Data.ThreeCycle, ThreeLength);
                Data.Ans3 += ThreeLength;
                memcpy(Data.Ans3, mpid3.str, mpid3.len);
                Data.Ans3 += mpid3.len;
                memcpy(Data.Ans3, mpid4.str, mpid4.len);
                Data.Ans3 += mpid4.len;
                memcpy(Data.Ans3, mpid5.str, mpid5.len);
                *(Data.Ans3 + mpid5.len - 1) = '\n';
                Data.Ans3 += mpid5.len;
                ++ans;
                continue;
              }
              const U32 &w7 = Data.LastWeight[v6];
              if (v1 == v6 || v2 == v6 || v3 == v6 || v4 == v6 ||
                  !judge(w6, w7) || !judge(w7, w1)) {
                continue;
              }
              const auto &mpid6 = MapID[v6];
              memcpy(Data.Ans4, Data.ThreeCycle, ThreeLength);
              Data.Ans4 += ThreeLength;
              memcpy(Data.Ans4, mpid3.str, mpid3.len);
              Data.Ans4 += mpid3.len;
              memcpy(Data.Ans4, mpid4.str, mpid4.len);
              Data.Ans4 += mpid4.len;
              memcpy(Data.Ans4, mpid5.str, mpid5.len);
              Data.Ans4 += mpid5.len;
              memcpy(Data.Ans4, mpid6.str, mpid6.len);
              *(Data.Ans4 + mpid6.len - 1) = '\n';
              Data.Ans4 += mpid6.len;
              ++ans;
            }
          }
        }
      }
    }
  }
  Data.answers += ans;
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
  Data.Ans0 = Data.Cycle0;
  Data.Ans1 = Data.Cycle1;
  Data.Ans2 = Data.Cycle2;
  Data.Ans3 = Data.Cycle3;
  Data.Ans4 = Data.Cycle4;
  while (true) {
    GetNextJob(job);
    if (job == -1) break;
    BackSearch(Data, job);
    ForwardSearch(Data, job);
  }
}

void SaveAnswer() {
  char firBuf[NUMLENGTH];
  U32 firIdx = NUMLENGTH;
  firBuf[--firIdx] = '\n';
  U32 x = Answers;
  if (x == 0) {
    firBuf[--firIdx] = '0';
  } else {
    while (x) {
      firBuf[--firIdx] = x % 10 + '0';
      x /= 10;
    }
  }
  FILE *fp = fopen(RESULT, "w");
  fwrite(firBuf + firIdx, 1, NUMLENGTH - firIdx, fp);
  fclose(fp);
}

void Simulation() {
  LoadData();
  std::thread Th[NTHREAD];
  for (int i = 0; i < NTHREAD; ++i) {
    Th[i] = std::thread(FindCircle, i);
  }
  for (auto &it : Th) it.join();
  for (auto &it : ThreadData) Answers += it.answers;
  std::cerr << "find over\n";
  // SaveAnswer();
#ifdef LOCAL
  std::cerr << "@ answers: " << Answers << "\n";
  for (int i = 0; i < NTHREAD; ++i) {
    std::cerr << "@ thread " << i << ": " << ThreadData[i].answers << "\n";
  }
#endif
}

int main() {
  Simulation();
  return 0;
}