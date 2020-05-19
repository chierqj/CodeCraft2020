#include <fcntl.h>
#include <sys/mman.h>
#include <sys/time.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <chrono>
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
#include <set>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#define uint uint32_t
#define ulong uint64_t
#define P10(x) ((x << 3) + (x << 1))

#ifdef LOCAL
#define TRAIN "../data/std/test_data.txt"
#define RESULT "../data/std/result.txt"
#else
#define TRAIN "/data/test_data.txt"
#define RESULT "/projects/student/result.txt"
#endif

namespace Color {
inline void reset() { fprintf(stderr, "\033[0m"); }
inline void red() { fprintf(stderr, "\033[1;31m"); }
inline void green() { fprintf(stderr, "\033[1;32m"); }
inline void yellow() { fprintf(stderr, "\033[1;33m"); }
inline void blue() { fprintf(stderr, "\033[1;34m"); }
inline void magenta() { fprintf(stderr, "\033[1;35m"); }
inline void cyan() { fprintf(stderr, "\033[1;36m"); }
inline void orange() { fprintf(stderr, "\033[38;5;214m"); }
inline void newline() { fprintf(stderr, "\n"); }
}  // namespace Color
class Timer {
 public:
  Timer() : m_begin(std::chrono::high_resolution_clock::now()) {}
  inline double elapsed() const {
    auto t = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - m_begin);
    double elapsed = (double)(t.count() * 1.0) / 1000.0;
    return elapsed;
  }

 private:
  std::chrono::time_point<std::chrono::high_resolution_clock> m_begin;
};

/*
 * 常量定义
 */
const uint MAX_EDGE = 3000000 + 7;  // 最多边数目
const uint MAX_NODE = 6000000 + 7;  // 最多点数目
const uint T = 8;                   // 线程个数

/************************** LoadData Begin ****************************/
struct HashTable {
  /*
   * 开源HashTable
   */
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
};
struct DFSEdge {
  uint idx;
  ulong w, val;  // 边权, 边宽; 目前默认边宽为1，等一手改需求
};
struct ReadEdge {
  uint u, v;
  ulong w, val;  // 边权, 边宽; 目前默认边宽为1，等一手改需求
};
struct LoadInfo {
  uint offsz = 0;                     // 线程hashmap偏移
  HashTable hashmap;                  // 线程独立hashmap
  std::vector<uint> ids;              // u % T = pid 的所有id
  std::vector<uint> th_ids[T];        // u % T的id
  std::vector<ReadEdge> edges;        // u % T = pid 的所有边
  std::vector<ReadEdge> th_edges[T];  // u % T 的边
};
uint g_NodeNum = 0;                      // 加载完数据后结点数目
uint g_EdgeNum = 0;                      // 加载完数据后边数目
uint Head[MAX_NODE], HeadLen[MAX_NODE];  // 前向星标记位置和长度
uint Back[MAX_NODE], BackLen[MAX_NODE];  // 后向星标记位置和长度
ReadEdge Edges[MAX_EDGE];                // 多线程边合并后的所有边
DFSEdge GHead[MAX_NODE];                 // 前向星所有边
DFSEdge GBack[MAX_NODE];                 // 后向星所有边
LoadInfo LoadInfos[T];                   // 多线程加载数据

/*
 * 加载数据
 * 1. 多线程解析buffer, 每个线程存 u%T 的边
 * 2. 将每个线程th_edges[u%T]的边memcpy到LoadInfo[pid]里面,并排序
 * 3. LoadInfo[pid]里面的边进行Hash
 * 4. 计算每个线程HashTable的偏移量
 * 5. 多线程对每个线程的HashTable重新Hash
 * 6. 合并HashTable
 * 7. 根据真实ID大小重新对HashTable的val进行排名, Rank
 * 8. 遍历所有边,将u, v改成Hash之后的映射id
 * 9. 构造前向星
 * 10. 构造后向星
 */
void InitLoadInfo() {
  for (uint i = 0; i < T; ++i) {
    LoadInfos[i].ids.reserve(MAX_NODE);
    LoadInfos[i].edges.reserve(MAX_EDGE);
    for (uint j = 0; j < T; ++j) {
      LoadInfos[i].th_ids[j].reserve(MAX_NODE);
      LoadInfos[i].th_edges[j].reserve(MAX_EDGE);
    }
  }
};
void addEdge(const uint &u, const uint &v, const ulong &w, const ulong &val,
             LoadInfo &data) {
  uint umod = u % T, vmod = v % T;
  data.th_edges[umod].emplace_back(ReadEdge{u, v, w, val});
  data.th_ids[umod].emplace_back(u);
  data.th_ids[vmod].emplace_back(v);
}
void HandleReadBuffer(const char *buffer, uint st, uint ed, uint pid) {
  const char *ptr = buffer + st, *end = buffer + ed;
  uint u = 0, v = 0;
  ulong w = 0, val = 1;
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
    while (*ptr != '\r' && *ptr != '\n') {
      w = P10(w) + *ptr - '0';
      ++ptr;
    }
    if (*ptr == '\r') ++ptr;
    ++ptr;
    addEdge(u, v, w, val, loadinfo);
    u = v = w = 0;
  }
}
void ReadBuffer() {
  uint fd = open(TRAIN, O_RDONLY);
  uint bufsize = lseek(fd, 0, SEEK_END);
  char *buffer = (char *)mmap(NULL, bufsize, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);
  std::thread Th[T];
  uint st = 0, block = bufsize / T;
  uint up = 0;
  for (uint i = 0; i < T; ++i) {
    if (i == T - 1) {
      Th[i] = std::thread(HandleReadBuffer, buffer, st, bufsize, i);
      up = T;
      break;
    }
    uint ed = st + block;
    while (buffer[ed] != '\n') ++ed;
    ++ed;
    Th[i] = std::thread(HandleReadBuffer, buffer, st, ed, i);
    st = ed;
    if (st >= bufsize) {
      up = i + 1;
      break;
    }
  }
  for (uint i = 0; i < up; ++i) Th[i].join();
}
void MergeSortIdAndEdge(uint pid) {
  auto &data = LoadInfos[pid];
  for (uint i = 0; i < T; ++i) {
    const auto &info = LoadInfos[i];
    if (!info.th_edges[pid].empty()) {
      data.edges.insert(data.edges.end(), info.th_edges[pid].begin(),
                        info.th_edges[pid].end());
    }
    if (!info.th_ids[pid].empty()) {
      data.ids.insert(data.ids.end(), info.th_ids[pid].begin(),
                      info.th_ids[pid].end());
    }
  }
  std::sort(data.ids.begin(), data.ids.end());
  std::sort(data.edges.begin(), data.edges.end(),
            [&](const ReadEdge &e1, const ReadEdge &e2) {
              if (e1.u == e2.u) return e1.v < e2.v;
              return e1.u < e2.u;
            });
}
void HashIds(uint pid) {
  auto &data = LoadInfos[pid];
  if (data.ids.empty()) return;
  uint pre = data.ids[0] + 7;
  for (auto &v : data.ids) {
    if (v == pre) continue;
    pre = v;
    data.hashmap.Insert(v);
  }
}
std::vector<uint> IDDom;
void CalMapOffSet() {
  for (uint i = 1; i < T; ++i) {
    const uint &sz = LoadInfos[i - 1].hashmap.Size();
    LoadInfos[i].offsz = LoadInfos[i - 1].offsz + sz;
  }
  g_NodeNum = LoadInfos[T - 1].offsz + LoadInfos[T - 1].hashmap.Size();
  IDDom.reserve(g_NodeNum);
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
  std::vector<uint> vec(g_NodeNum);
  Rank.reserve(g_NodeNum);
  for (uint i = 0; i < g_NodeNum; ++i) vec[i] = i;
  std::sort(vec.begin(), vec.end(),
            [&](const uint &x, const uint &y) { return IDDom[x] < IDDom[y]; });
  for (uint i = 0; i < g_NodeNum; ++i) {
    const uint &x = vec[i];
    Rank[x] = i;
  }
}
void ReHashTable(uint pid) {
  auto &hashmap = LoadInfos[pid].hashmap;
  for (uint i = 0; i < hashmap.Size(); ++i) {
    auto &p = hashmap.Map[hashmap.HashIdx[i]];
    p.val = Rank[p.val];
    IDDom[p.val] = p.key;
  }
}
void RankEdge(uint pid) {
  auto &data = LoadInfos[pid];
  for (auto &e : data.edges) {
    e.u = LoadInfos[e.u % T].hashmap.Query(e.u);
    e.v = LoadInfos[e.v % T].hashmap.Query(e.v);
  }
}
void MergeEdge() {
  uint left[T] = {0};
  ReadEdge *ptr = Edges;
  uint cnt = 0;
  while (true) {
    uint minx = g_NodeNum + 7, minidx = g_NodeNum + 7;
    for (uint i = 0; i < T; ++i) {
      const auto &data = LoadInfos[i];
      if (left[i] < data.edges.size() && data.edges[left[i]].u < minx) {
        minx = data.edges[left[i]].u;
        minidx = i;
      }
    }
    if (minidx == g_NodeNum + 7) break;
    const auto &data = LoadInfos[minidx];
    uint &l = left[minidx], r = data.edges.size();
    while (l < r && data.edges[l].u == minx) {
      memcpy(ptr, &data.edges[left[minidx]], sizeof(ReadEdge));
      ++ptr;
      ++l;
      ++cnt;
    }
  }
  g_EdgeNum = ptr - Edges;
}
void BuildGraphHead(uint pid) {
  uint pre = g_NodeNum + 7;
  for (uint i = 0; i < g_EdgeNum; ++i) {
    const auto &e = Edges[i];
    if (e.u % T == pid) {
      if (e.u != pre) {
        Head[e.u] = i;
      }
      ++HeadLen[e.u];
      pre = e.u;
      GHead[i] = {e.v, e.w, e.val};
    }
    if (e.v % T == pid) {
      ++BackLen[e.v];
    }
  }
}
void BuildGraphBack(uint pid) {
  std::vector<uint> cnt(g_NodeNum, 0);
  for (int i = g_EdgeNum - 1; i >= 0; --i) {
    const auto &e = Edges[i];
    if (e.v % T == pid) {
      auto &p = GBack[Back[e.v] + cnt[e.v]];
      p = {e.u, e.w, e.val};
      ++cnt[e.v];
    }
  }
}
void LoadData() {
#ifdef LOCAL
  Timer t;
#endif
  InitLoadInfo();
  ReadBuffer();
  std::thread Th[T];
  for (uint i = 0; i < T; ++i) Th[i] = std::thread(MergeSortIdAndEdge, i);
  for (auto &it : Th) it.join();
  for (uint i = 0; i < T; ++i) Th[i] = std::thread(HashIds, i);
  for (auto &it : Th) it.join();
  CalMapOffSet();
  for (uint i = 0; i < T; ++i) Th[i] = std::thread(MergeHashTable, i);
  for (auto &it : Th) it.join();
  RankHashTable();
  for (uint i = 0; i < T; ++i) Th[i] = std::thread(ReHashTable, i);
  for (auto &it : Th) it.join();
  for (uint i = 0; i < T; ++i) Th[i] = std::thread(RankEdge, i);
  for (auto &it : Th) it.join();
  MergeEdge();
  for (uint i = 0; i < T; ++i) Th[i] = std::thread(BuildGraphHead, i);
  for (auto &it : Th) it.join();
  for (uint i = 1; i <= g_NodeNum; ++i) {
    Head[i] = Head[i - 1] + HeadLen[i - 1];
    Back[i] = Back[i - 1] + BackLen[i - 1];
  }
  for (uint i = 0; i < T; ++i) Th[i] = std::thread(BuildGraphBack, i);
  for (auto &it : Th) it.join();
#ifdef LOCAL
  Color::green();
  printf("@ LoadData: [结点: %d, 边: %d] [cost: %.4fs]\n", g_NodeNum, g_EdgeNum,
         t.elapsed());
  Color::reset();
#endif
}
/************************** LoadData End ****************************/

/************************** Algorithm Start ****************************/
typedef std::pair<uint, double> prud;
struct Node {
  uint idx;   // 结点id
  ulong val;  // dis[idx]
  bool operator<(const Node &r) const { return val > r.val; }
};
struct ThreadData {
  std::priority_queue<Node> pq;  // 优先队列
  std::queue<uint> q;            // SPFA
  std::vector<bool> vis;         // 标记访问
  std::vector<ulong> dis;        // 最短距离
  std::vector<uint> count;       // (dis[u] + w == dis[v]) -> count[v]
  uint points[MAX_NODE];         // 最短路上的点
  std::vector<double> ans;       // 线程结果
};
ThreadData TData[T];       // 线程数据
std::vector<prud> Answer;  // 最终结果

/*
 * 1. dijkstr寻找最短路
 * 2. 利用拆分的公式计算中心性
 */
void ClearThreadData(ThreadData &Data) {
  Data.vis = std::vector<bool>(g_NodeNum);
  Data.dis = std::vector<ulong>(g_NodeNum, std::numeric_limits<ulong>::max());
  Data.count = std::vector<uint>(g_NodeNum, 0);
}

void SPFAAnswer(ThreadData &Data, uint start) {
  uint l = 0, r = 0;
  std::vector<double> f(g_NodeNum, 0);
  Data.points[++r] = start;
  f[start] = 1;
  // 先计算拓扑序
  while (r > l) {
    int u = Data.points[++l];
    for (uint i = Head[u]; i < Head[u + 1]; ++i) {
      const auto &e = GHead[i];
      if (Data.dis[u] + e.w == Data.dis[e.idx]) {
        if (!--Data.count[e.idx]) Data.points[++r] = e.idx;
        f[e.idx] += f[u];
      }
    }
  }
  // 更新ans
  std::vector<double> g(g_NodeNum, 0);
  for (int p = r; p > 1; --p) {
    uint u = Data.points[p];
    for (uint i = Head[u]; i < Head[u + 1]; ++i) {
      const auto &e = GHead[i];
      if (Data.dis[u] + e.w == Data.dis[e.idx]) {
        g[u] += g[e.idx];
      }
    }
    Data.ans[u] += 1.0 * f[u] * g[u];
    g[u] += (double)(1.0 / f[u]);
  }
}
void SPFA(ThreadData &Data, uint start) {
  Data.q.push(start);
  Data.vis[start] = true;
  Data.dis[start] = 0;
  while (!Data.q.empty()) {
    const uint u = Data.q.front();
    Data.q.pop();
    Data.vis[u] = false;
    for (uint i = Head[u]; i < Head[u + 1]; ++i) {
      const auto &e = GHead[i];
      if (Data.dis[u] + e.w > Data.dis[e.idx]) continue;
      if (Data.dis[u] + e.w == Data.dis[e.idx]) {
        ++Data.count[e.idx];  // (dis[u] + w = dis[v]) -> count[v]++;
      } else {
        Data.dis[e.idx] = Data.dis[u] + e.w;
        Data.count[e.idx] = 1;
        if (!Data.vis[e.idx]) {
          Data.q.push(e.idx);
          Data.vis[e.idx] = true;
        }
      }
    }
  }
  SPFAAnswer(Data, start);
}
void Dijkstra(ThreadData &Data, uint start) {
  Data.pq.push(Node{start, 0});
  uint l = 0, r = 0;
  Data.dis[start] = 0;
  Data.count[start] = 1;
  while (!Data.pq.empty()) {
    const auto u = Data.pq.top().idx;
    Data.pq.pop();
    if (Data.vis[u]) continue;
    Data.vis[u] = true;
    Data.points[++r] = u;  //拓扑序
    for (uint i = Head[u]; i < Head[u + 1]; ++i) {
      const auto &e = GHead[i];
      if (Data.dis[u] + e.w > Data.dis[e.idx]) continue;
      if (Data.dis[u] + e.w == Data.dis[e.idx]) {
        Data.count[e.idx] += Data.count[u];  // s到e.idx的最短路条数
      } else {
        Data.count[e.idx] = Data.count[u];
        Data.dis[e.idx] = Data.dis[u] + e.w;
        Data.pq.push(Node{e.idx, Data.dis[e.idx]});
      }
    }
  }

  //更新ans
  std::vector<double> g(g_NodeNum, 0);
  for (int p = r; p > 1; --p) {
    uint u = Data.points[p];
    for (uint i = Head[u]; i < Head[u + 1]; ++i) {
      const auto &e = GHead[i];
      if (Data.dis[u] + e.w == Data.dis[e.idx]) {
        g[u] += g[e.idx];
      }
    }
    Data.ans[u] += 1.0 * Data.count[u] * g[u];
    g[u] += (double)(1.0 / Data.count[u]);
  }
}

// atomic for job
inline void getJob(uint &job) {
  static uint l = 0;
  static std::atomic_flag lock = ATOMIC_FLAG_INIT;
  while (lock.test_and_set())
    ;
  job = l < g_NodeNum ? l++ : -1;
  lock.clear();
}

// 1: dijkstra; 2: spfa
inline uint ChooseStragety(const uint &start) { return random() & 1; }

void HandleSolve(uint pid) {
  auto &Data = TData[pid];
  Data.ans = std::vector<double>(g_NodeNum, 0);

  ClearThreadData(Data);
  uint job = 0;
  while (true) {
    getJob(job);
    if (job == -1) break;

    uint stragety = ChooseStragety(job);

    if (stragety == 1) {
      Dijkstra(Data, job);
    } else {
      SPFA(Data, job);
    }
    ClearThreadData(Data);
  }
}

void Solve() {
#ifdef LOCAL
  Timer t;
#endif

  std::thread Th[T];
  for (uint i = 0; i < T; ++i) Th[i] = std::thread(HandleSolve, i);
  for (auto &it : Th) it.join();

  Answer = std::vector<prud>(g_NodeNum);
  for (uint i = 0; i < g_NodeNum; ++i) {
    Answer[i].first = IDDom[i];
    Answer[i].second = 0;
    for (const auto &data : TData) {
      Answer[i].second += data.ans[i];
    }
  }

#ifdef LOCAL
  Color::magenta();
  printf("@ Solve: [dijkstra] [cost: %.4fs]\n", t.elapsed());
  Color::reset();
#endif
}

/************************** Algorithm End ****************************/

void SaveAnswer() {
#ifdef LOCAL
  Timer t;
#endif
  std::sort(Answer.begin(), Answer.end(), [](const prud &p1, const prud &p2) {
    if (p1.second == p2.second) return p1.first < p2.first;
    return p1.second > p2.second;
  });

  uint up = g_NodeNum > 100 ? 100 : g_NodeNum;

  char *buffer = new char[50 * 100];
  char *ptr = buffer;

  for (uint i = 0; i < up; ++i) {
    const auto &it = Answer[i];
    ptr += sprintf(ptr, "%d,%.3lf\n", it.first, it.second);
  }

  FILE *fp = fopen(RESULT, "w");
  fwrite(buffer, 1, ptr - buffer, fp);
  fclose(fp);

#ifdef LOCAL
  Color::magenta();
  printf("@ Save: [sort and save] [cost: %.4fs]\n", t.elapsed());
  Color::reset();
#endif
}

int main() {
  std::cerr << std::fixed << std::setprecision(3);

  LoadData();
  Solve();
  SaveAnswer();

  return 0;
}