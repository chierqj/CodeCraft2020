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
const uint MAXEDGE = 3000000 + 7;    // 最多边数目
const uint MAXN = MAXEDGE << 1;      // 最多点数目
const uint MAXCYCLE = 20000000 + 7;  // 最多环数
const uint NTHREAD = 4;              // 线程个数
const uint NUMLENGTH = 11;           // ID最大长度
uint ST[NTHREAD], ED[NTHREAD];

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
  uint ReachPointCount = 0;  // 反向可达点数目
  char Reach[MAXN];          // 标记反向可达
  uint ReachPoint[MAXN];     // 可达点集合
  uint LastWeight[MAXN];     // 最后一步权重
  char *Cycle0 = new char[MAXCYCLE * NUMLENGTH * 3];
  char *Cycle1 = new char[MAXCYCLE * NUMLENGTH * 4];
  char *Cycle2 = new char[MAXCYCLE * NUMLENGTH * 5];
  char *Cycle3 = new char[MAXCYCLE * NUMLENGTH * 6];
  char *Cycle4 = new char[MAXCYCLE * NUMLENGTH * 7];
  char *ans0, *ans1, *ans2, *ans3, *ans4;
  uint bufsize0, bufsize1, bufsize2, bufsize3, bufsize4;
  uint bufsize;
} ThreadData[NTHREAD];  // 线程找环

/*
 * 答案
 */
char FirstBuf[NUMLENGTH];
uint FirstBufLen;
uint TotalBufferSize = 0;

/*
 * atomic 锁
 */
uint EdgeCur = 0;
uint IDCur = 0;
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
  uint ans = 0;
  const auto &children1 = Children[st];
  const auto &mpid0 = MapID[st];
  for (const auto &it1 : children1) {
    const uint &v1 = it1.first, &w1 = it1.second;
    if (v1 < st) continue;
    const auto &children2 = Children[v1];
    const auto &mpid1 = MapID[v1];
    for (const auto &it2 : children2) {
      const uint &v2 = it2.first, &w2 = it2.second;
      if (v2 <= st || !judge(w1, w2)) continue;
      const auto &children3 = Children[v2];
      const auto &mpid2 = MapID[v2];
      for (const auto &it3 : children3) {
        const uint &v3 = it3.first, &w3 = it3.second;
        if (v3 < st || v3 == v1 || !judge(w2, w3)) {
          continue;
        } else if (v3 == st) {
          if (!judge(w3, w1)) continue;
          memcpy(Data.ans0, mpid0.str, mpid0.len);
          Data.ans0 += mpid0.len;
          memcpy(Data.ans0, mpid1.str, mpid1.len);
          Data.ans0 += mpid1.len;
          memcpy(Data.ans0, mpid2.str, mpid2.len);
          Data.ans0[mpid2.len - 1] = '\n';
          Data.ans0 += mpid2.len;
          ++ans;
          continue;
        }
        const auto &children4 = Children[v3];
        const auto &mpid3 = MapID[v3];
        for (const auto &it4 : children4) {
          const uint &v4 = it4.first, &w4 = it4.second;
          if (!(Data.Reach[v4] & 4) || !judge(w3, w4)) {
            continue;
          } else if (v4 == st) {
            if (!judge(w4, w1)) continue;
            memcpy(Data.ans1, mpid0.str, mpid0.len);
            Data.ans1 += mpid0.len;
            memcpy(Data.ans1, mpid1.str, mpid1.len);
            Data.ans1 += mpid1.len;
            memcpy(Data.ans1, mpid2.str, mpid2.len);
            Data.ans1 += mpid2.len;
            memcpy(Data.ans1, mpid3.str, mpid3.len);
            Data.ans1[mpid3.len - 1] = '\n';
            Data.ans1 += mpid3.len;
            ++ans;
            continue;
          } else if (v1 == v4 || v2 == v4) {
            continue;
          }
          const auto &children5 = Children[v4];
          const auto &mpid4 = MapID[v4];
          for (const auto &it5 : children5) {
            const uint &v5 = it5.first, &w5 = it5.second;
            if (!Data.Reach[v5] || !judge(w4, w5)) {
              continue;
            } else if (v5 == st) {
              if (!judge(w5, w1)) continue;
              memcpy(Data.ans2, mpid0.str, mpid0.len);
              Data.ans2 += mpid0.len;
              memcpy(Data.ans2, mpid1.str, mpid1.len);
              Data.ans2 += mpid1.len;
              memcpy(Data.ans2, mpid2.str, mpid2.len);
              Data.ans2 += mpid2.len;
              memcpy(Data.ans2, mpid3.str, mpid3.len);
              Data.ans2 += mpid3.len;
              memcpy(Data.ans2, mpid4.str, mpid4.len);
              Data.ans2[mpid4.len - 1] = '\n';
              Data.ans2 += mpid4.len;
              ++ans;
              continue;
            } else if (v1 == v5 || v2 == v5 || v3 == v5) {
              continue;
            }
            const auto &children6 = Children[v5];
            const auto &mpid5 = MapID[v5];
            for (const auto &it6 : children6) {
              const uint &v6 = it6.first, &w6 = it6.second;
              if (!(Data.Reach[v6] & 1) || !judge(w5, w6)) {
                continue;
              } else if (v6 == st) {
                if (!judge(w6, w1)) continue;
                memcpy(Data.ans3, mpid0.str, mpid0.len);
                Data.ans3 += mpid0.len;
                memcpy(Data.ans3, mpid1.str, mpid1.len);
                Data.ans3 += mpid1.len;
                memcpy(Data.ans3, mpid2.str, mpid2.len);
                Data.ans3 += mpid2.len;
                memcpy(Data.ans3, mpid3.str, mpid3.len);
                Data.ans3 += mpid3.len;
                memcpy(Data.ans3, mpid4.str, mpid4.len);
                Data.ans3 += mpid4.len;
                memcpy(Data.ans3, mpid5.str, mpid5.len);
                Data.ans3[mpid5.len - 1] = '\n';
                Data.ans3 += mpid5.len;

                ++ans;
                continue;
              }
              const uint &w7 = Data.LastWeight[v6];
              if (v1 == v6 || v2 == v6 || v3 == v6 || v4 == v6 ||
                  !judge(w6, w7) || !judge(w7, w1)) {
                continue;
              }
              const auto &mpid6 = MapID[v6];
              memcpy(Data.ans4, mpid0.str, mpid0.len);
              Data.ans4 += mpid0.len;
              memcpy(Data.ans4, mpid1.str, mpid1.len);
              Data.ans4 += mpid1.len;
              memcpy(Data.ans4, mpid2.str, mpid2.len);
              Data.ans4 += mpid2.len;
              memcpy(Data.ans4, mpid3.str, mpid3.len);
              Data.ans4 += mpid3.len;
              memcpy(Data.ans4, mpid4.str, mpid4.len);
              Data.ans4 += mpid4.len;
              memcpy(Data.ans4, mpid5.str, mpid5.len);
              Data.ans4 += mpid5.len;
              memcpy(Data.ans4, mpid6.str, mpid6.len);
              Data.ans4[mpid6.len - 1] = '\n';
              Data.ans4 += mpid6.len;
              ++ans;
            }
          }
        }
      }
    }
  }
  Data.answers += ans;
}

void FindCircle(uint pid) {
  auto &Data = ThreadData[pid];
  Data.ans0 = Data.Cycle0;
  Data.ans1 = Data.Cycle1;
  Data.ans2 = Data.Cycle2;
  Data.ans3 = Data.Cycle3;
  Data.ans4 = Data.Cycle4;
  uint st = ST[pid], ed = ED[pid];
  for (uint i = st; i < ed; ++i) {
    const uint &job = Jobs[i];
    BackSearch(Data, job);
    ForwardSearch(Data, job);
  }
  Data.bufsize0 = Data.ans0 - Data.Cycle0;
  Data.bufsize1 = Data.ans1 - Data.Cycle1;
  Data.bufsize2 = Data.ans2 - Data.Cycle2;
  Data.bufsize3 = Data.ans3 - Data.Cycle3;
  Data.bufsize4 = Data.ans4 - Data.Cycle4;
  Data.bufsize = Data.bufsize0 + Data.bufsize1 + Data.bufsize2 + Data.bufsize3 +
                 Data.bufsize4;
  TotalBufferSize += Data.bufsize;
}

uint OffSet[NTHREAD][5];
void CalOffSet() {
  uint x = FirstBufLen;
  for (uint i = 0; i < NTHREAD; ++i) {
    OffSet[i][0] = x;
    x += ThreadData[i].bufsize0;
  }
  for (uint i = 0; i < NTHREAD; ++i) {
    OffSet[i][1] = x;
    x += ThreadData[i].bufsize1;
  }
  for (uint i = 0; i < NTHREAD; ++i) {
    OffSet[i][2] = x;
    x += ThreadData[i].bufsize2;
  }
  for (uint i = 0; i < NTHREAD; ++i) {
    OffSet[i][3] = x;
    x += ThreadData[i].bufsize3;
  }
  for (uint i = 0; i < NTHREAD; ++i) {
    OffSet[i][4] = x;
    x += ThreadData[i].bufsize4;
  }
#ifdef TEST
  for (int i = 0; i < NTHREAD; ++i) {
    for (int j = 0; j < 5; ++j) {
      std::cerr << OffSet[i][j] << ", ";
    }
    std::cerr << "\n";
  }
#endif
}
void HandleSaveAnswer(int pid) {
  uint fd = open(RESULT, O_RDWR | O_CREAT, 0666);
  char *result = (char *)mmap(NULL, TotalBufferSize, PROT_READ | PROT_WRITE,
                              MAP_SHARED, fd, 0);
  char *ptr = result;
  ftruncate(fd, TotalBufferSize);
  close(fd);
  if (pid == 0) {
    memcpy(result, FirstBuf, FirstBufLen);
  }
  const auto &data = ThreadData[pid];
  result = ptr + OffSet[pid][0];
  memcpy(result, data.Cycle0, data.bufsize0);
  result = ptr + OffSet[pid][1];
  memcpy(result, data.Cycle1, data.bufsize1);
  result = ptr + OffSet[pid][2];
  memcpy(result, data.Cycle2, data.bufsize2);
  result = ptr + OffSet[pid][3];
  memcpy(result, data.Cycle3, data.bufsize3);
  result = ptr + OffSet[pid][4];
  memcpy(result, data.Cycle4, data.bufsize4);
}
void SaveAnswer() {
  sprintf(FirstBuf, "%d\n", Answers);
  FirstBufLen = strlen(FirstBuf);
  TotalBufferSize += FirstBufLen;
  CalOffSet();
  std::thread Th[NTHREAD];
  for (int i = 1; i < NTHREAD; ++i) Th[i] = std::thread(HandleSaveAnswer, i);
  HandleSaveAnswer(0);
  for (int i = 1; i < NTHREAD; ++i) Th[i].join();
}

void SetParam() {
  double sum = 0;
  double DC = JobsCount;
  double pre[JobsCount];
  for (int i = 0; i < JobsCount; ++i) {
    sum += std::pow((1.0 - i / DC), 3.6) * Children[Jobs[i]].size();
    pre[i] = sum;
  }
  double block = sum / NTHREAD;
  double val = block;
  for (int i = 0; i < NTHREAD - 1; ++i) {
    int p = std::lower_bound(pre, pre + JobsCount, val) - pre;
    val += block;
    ED[i] = p;
    ST[i + 1] = p;
  }
  ED[NTHREAD - 1] = JobsCount;
#ifdef LOCAL
  for (int i = 0; i < NTHREAD; ++i) {
    std::cerr << "@ param: (" << ST[i] << ", " << ED[i] << ")\n";
  }
#endif
}
void Simulation() {
  LoadData();
  SetParam();
  std::thread Th[NTHREAD];
  for (uint i = 1; i < NTHREAD; ++i) Th[i] = std::thread(FindCircle, i);
  FindCircle(0);
  for (uint i = 1; i < NTHREAD; ++i) Th[i].join();
  for (auto &it : ThreadData) Answers += it.answers;
#ifdef LOCAL
  uint idx = 0;
  for (auto &it : ThreadData) {
    std::cerr << "@ thread " << idx++ << ": " << it.answers << ", "
              << it.bufsize << "\n";
  }
  std::cerr << "@ answers: " << Answers << ", " << TotalBufferSize << "\n";
#endif
  SaveAnswer();
}

int main() {
  Simulation();
  return 0;
}