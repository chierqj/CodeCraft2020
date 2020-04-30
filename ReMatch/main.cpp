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
#define U32 uint32_t

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
} __attribute__((packed)) MapID[MAXNODE];
struct ThData {
  U32 answers;
  char Reachable[MAXNODE];
  std::vector<U32> ReachablePoint;
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
std::vector<U32> Parents[MAXNODE];
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
    E[count++] = std::make_tuple(p2, p2, ew);
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
    while (*ptr != '\n') {
      w = w * 10 + *ptr - '0';
      ++ptr;
    }
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
    for (int j = 0; j < ThEdgesCount[i]; ++j) {
      const auto &e = ThEdges[i][j];
      const U32 &eu = std::get<0>(e);
      const U32 &ev = std::get<1>(e);
      const U32 &ew = std::get<2>(e);
      Children[eu].emplace_back(std::make_pair(ev, ew));
      Parents[ev].emplace_back(eu);
    }
  }
  for (int i = 0; i < MaxID; ++i) {
    if (!Children[i].empty() && !Children[i].empty()) {
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
  for (auto &v : Data.ReachablePoint) Data.Reachable[v] = 0;
  Data.ReachablePoint.clear();
  Data.ReachablePoint.emplace_back(job);
  Data.Reachable[job] = 7;
  for (auto &it1 : Parents[job]) {
    if (it1 <= job) continue;
    Data.Reachable[it1] = 7;
    Data.ReachablePoint.emplace_back(it1);
    for (auto &it2 : Parents[it1]) {
      if (it2 <= job) continue;
      Data.Reachable[it2] |= 6;
      Data.ReachablePoint.emplace_back(it2);
      for (auto &it3 : Parents[it2]) {
        if (it3 <= job || it3 == it1) continue;
        Data.Reachable[it3] |= 4;
        Data.ReachablePoint.emplace_back(it3);
      }
    }
  }
}
/*
void ForwardSearch(ThData &Data, const U32 &job) {
  auto &ans = Answer[job];
  for (auto &it1 : Children[job]) {
    if (it1.first < job) continue;
    for (auto &it2 : Children[it1.first]) {
      if (it2.first <= job) continue;
      for (auto &it3 : Children[it2.first]) {
        if (it3.first < job || it3.first == it1.first) continue;
        if (it3.first == job) {
          ++Data.answers;
          ans[0].emplace_back(it1.first);
          ans[0].emplace_back(it2.first);
          continue;
        }
        for (auto &it4 : Children[it3.first]) {
          if (!(Data.Reachable[it4.first] & 4)) continue;
          if (it4.first == job) {
            ++Data.answers;
            ans[1].emplace_back(it1.first);
            ans[1].emplace_back(it2.first);
            ans[1].emplace_back(it3.first);
            continue;
          }
          if (it4.first == it1.first || it4.first == it2.first) {
            continue;
          }
          for (auto &it5 : Children[it4.first]) {
            if (!(Data.Reachable[it5.first] & 2)) continue;
            if (it5.first == job) {
              ++Data.answers;
              ans[2].emplace_back(it1.first);
              ans[2].emplace_back(it2.first);
              ans[2].emplace_back(it3.first);
              ans[2].emplace_back(it4.first);
              continue;
            }
            if (it5.first == it1.first || it5.first == it2.first ||
                it5.first == it3.first) {
              continue;
            }
            for (auto &it6 : Children[it5.first]) {
              if (!(Data.Reachable[it6.first] & 1)) continue;
              if (it6.first == job) {
                ++Data.answers;
                ans[3].emplace_back(it1.first);
                ans[3].emplace_back(it2.first);
                ans[3].emplace_back(it3.first);
                ans[3].emplace_back(it4.first);
                ans[3].emplace_back(it5.first);
                continue;
              }
              if (it6.first == it1.first || it6.first == it2.first ||
                  it6.first == it3.first || it6.first == it4.first) {
                continue;
              }
              ++Data.answers;
              ans[4].emplace_back(it1.first);
              ans[4].emplace_back(it2.first);
              ans[4].emplace_back(it3.first);
              ans[4].emplace_back(it4.first);
              ans[4].emplace_back(it5.first);
              ans[4].emplace_back(it6.first);
            }
          }
        }
      }
    }
  }
}
*/

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

/*
void FindCircle(int pid) {
  U32 job = 0;
  auto &Data = ThreadData[pid];
  while (true) {
    GetNextJob(job);
    if (job == -1) break;
    Answer[job].resize(5);
    BackSearch(Data, job);
    ForwardSearch(Data, job);
  }
}
*/

/*
void SaveAnswer() {
  FILE *fp = fopen(RESULT, "w");
  char tmp[10];
  int tidx = 10;
  U32 tol = answers;
  tmp[--tidx] = '\n';
  while (tol) {
    tmp[--tidx] = tol % 10 + '0';
    tol /= 10;
  }
  fwrite(tmp + tidx, 1, 10 - tidx, fp);
  char line[1];
  line[0] = '\n';
  for (int sz = 0; sz < 5; ++sz) {
    for (auto &job : Jobs) {
      if (Answer[job].empty()) continue;
      const auto &ans = Answer[job][sz];
      if (ans.empty()) continue;
      const auto &mpjob = MapID[job];
      int idx = 0;
      for (auto &v : ans) {
        if (!idx) fwrite(mpjob.str, 1, mpjob.len, fp);
        ++idx;
        if (idx == sz + 2) {
          idx = 0;
          fwrite(MapID[v].str, 1, MapID[v].len - 1, fp);
          fwrite(line, 1, 1, fp);
        } else {
          fwrite(MapID[v].str, 1, MapID[v].len, fp);
        }
      }
    }
  }
  fclose(fp);
}
*/
void Simulation() {
  LoadData();
  /*
  std::thread Th[NTHREAD];
  for (int i = 0; i < NTHREAD; ++i) {
    Th[i] = std::thread(&FindCircle, this, i);
  }
  for (auto &it : Th) it.join();
  for (auto &it : ThreadData) answers += it.answers;
  SaveAnswer();
  */

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