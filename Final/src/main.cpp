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
#include <deque>
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
#define TRAIN "../data/std3/test_data.txt"
#define RESULT "../data/std3/result.txt"
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
const uint MAX_EDGE = 2500000 + 7;  // 最多边数目
const uint MAX_NODE = 5000000 + 7;  // 最多点数目
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
  uint m_count = 0;
  Data Map[MOD1];
  uint HashIdx[MOD1];
  uint Size() { return m_count; }
  uint hash(const uint &k, const int &i) {
    return (k % MOD1 + i * (MOD2 - k % MOD2)) % MOD1;
  }
  void Insert(const uint &key) {
    for (int i = 0; i < MOD1; ++i) {
      const uint &val = this->hash(key, i);
      if (Map[val].val == -1) {
        Map[val].key = key;
        Map[val].val = m_count;
        HashIdx[m_count++] = val;
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
  uint w;
};
struct ReadEdge {
  uint u, v;
  uint w;
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
 * 是否是稀疏图
 * 稀疏图: std::priority_queue
 * 稠密图: 手写heap
 */
bool IfSparseGraph = false;

/*
 * Team: 孤芳自赏
 * No1. chier
 * No2. XDUls
 * No3. yangzhi__
 *
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
 *
 * 注释不够多，不美观
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
void addEdge(const uint &u, const uint &v, const uint &w, LoadInfo &data) {
  if (w == 0) return;
  uint umod = u % T, vmod = v % T;
  data.th_edges[umod].emplace_back(ReadEdge{u, v, w});
  data.th_ids[umod].emplace_back(u);
  data.th_ids[vmod].emplace_back(v);
}
void HandleReadBuffer(const char *buffer, uint st, uint ed, uint pid) {
  const char *ptr = buffer + st, *end = buffer + ed;
  uint u = 0, v = 0;
  uint w = 0;
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
    addEdge(u, v, w, loadinfo);
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
  uint m_pre = data.ids[0] + 7;
  for (auto &v : data.ids) {
    if (v == m_pre) continue;
    m_pre = v;
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
  uint m_pre = g_NodeNum + 7;
  for (uint i = 0; i < g_EdgeNum; ++i) {
    const auto &e = Edges[i];
    if (e.u % T == pid) {
      if (e.u != m_pre) {
        Head[e.u] = i;
      }
      ++HeadLen[e.u];
      m_pre = e.u;
      GHead[i] = {e.v, e.w};
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
      p = {e.u, e.w};
      ++cnt[e.v];
    }
  }
}
void LoadData() {
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
}
/************************** LoadData End ****************************/

/************************** Algorithm Start ****************************/

/*
 * Team: 孤芳自赏
 * No1. chier
 * No2. XDUls
 * No3. yangzhi__
 *
 * [这里是DataStruct的定义]
 * 没啥说的，就是变量定义，都加注释了，自己看
 *
 * 注释不够多，不美观
 */

uint Fa[MAX_NODE];                       // 并查集
uint BlockNum = 0;                       // 集合数目
std::pair<uint, uint> Blocks[MAX_NODE];  // 并查集(V,E)
ulong Label[MAX_NODE];                   // topsort 标记
bool Top[MAX_NODE];                      // Top

uint getFa(uint x) { return Fa[x] == x ? x : Fa[x] = getFa(Fa[x]); }

struct SolverData {
  uint pointNum = 0;        // 拓扑点数目
  uint points[MAX_NODE];    // 拓扑点
  uint count[MAX_NODE];     // 最短路径数目
  uint viscount[MAX_NODE];  // viscount
  uint head[MAX_NODE];      // for heap
  uint id[MAX_NODE];        // for heap
  ulong dis[MAX_NODE];      // 最短距离
  bool vis[MAX_NODE];       // vis
  double ans[MAX_NODE];     // 保存答案
  double g[MAX_NODE];       // gvalue
};
SolverData SovData[T];

struct Node {
  uint idx;
  ulong val;
  bool operator<(const Node &r) const { return val > r.val; }
};

/*
 * 全地球人公用一个更新答案的接口
 * 需要正确计算的东西 (拓扑序点 + m_count)
 * 1. 按照公式更新。
 * 2. 如果是关键点，就反向bfs更新
 */
void GetAnswer(SolverData &Data, const uint &start) {
  uint(&m_points)[MAX_NODE] = Data.points;
  uint(&m_count)[MAX_NODE] = Data.count;
  ulong(&m_dis)[MAX_NODE] = Data.dis;
  double(&m_ans)[MAX_NODE] = Data.ans;
  double(&m_g)[MAX_NODE] = Data.g;
  uint &m_pointNum = Data.pointNum;
  // 更新ans
  double pw = Label[start] + 1;
  for (uint p = m_pointNum; p > 1; --p) {
    const uint &u = m_points[p];
    const DFSEdge *e = &GHead[Head[u]];
    const auto &l = Head[u], &r = Head[u + 1];
    for (uint i = l; i < r; ++i, ++e) {
      const uint &v = e->idx;
      if (m_dis[u] + e->w == m_dis[v]) {
        m_g[u] += m_g[v];
      }
    }
    m_ans[u] += (double)(m_g[u] * (double)m_count[u] * pw);
    m_g[u] += (double)(1.0 / (double)(m_count[u]));
  }

  if (Label[start] > 0) {
    m_ans[start] += (double)(m_pointNum - 1) * Label[start];
    std::queue<std::pair<uint, ulong>> q;
    q.push(std::make_pair(start, m_pointNum - 1));
    while (!q.empty()) {
      auto head = q.front();
      const uint &u = head.first;
      const ulong &cnt = head.second;
      q.pop();
      for (uint i = Back[u]; i < Back[u + 1]; ++i) {
        const auto &e = GBack[i];
        if (Label[e.idx] <= 0 || !Top[e.idx] || HeadLen[e.idx] != 1) continue;
        double x = Label[e.idx] * (cnt + 1);
        m_ans[e.idx] += x;
        q.push(std::make_pair(e.idx, cnt + 1));
      }
    }
  }
}

// 全地球人公用一个clear的接口
void Clear(SolverData &Data) {
  uint(&m_points)[MAX_NODE] = Data.points;
  uint(&m_count)[MAX_NODE] = Data.count;
  uint(&m_viscount)[MAX_NODE] = Data.viscount;
  uint(&m_id)[MAX_NODE] = Data.id;
  bool(&m_vis)[MAX_NODE] = Data.vis;
  ulong(&m_dis)[MAX_NODE] = Data.dis;
  double(&m_g)[MAX_NODE] = Data.g;
  uint &m_pointNum = Data.pointNum;
  for (uint i = 1; i <= m_pointNum; ++i) {
    const uint &v = m_points[i];
    m_dis[v] = UINT64_MAX;
    m_g[v] = 0;
    m_vis[v] = false;
    m_viscount[v] = 0;
    m_count[v] = 0;
    m_id[v] = 0;
  }
  m_pointNum = 0;
}

/*
 * Team: 孤芳自赏
 * No1. chier
 * No2. XDUls
 * No3. yangzhi__
 *
 * [Stragety: Dijkstra]
 *
 * 这里是一个Dijkstra算法，为了看起来好看，我必须加这个注释
 * 1. std::priority_queue 优先队列记录count
 * 2. 全地球人公用一个更新答案的接口
 * 3. 全地球人公用一个clear的接口
 *
 * 注释不够多，不美观
 */

void Dijkstra(SolverData &Data, const uint &start) {
  uint(&m_points)[MAX_NODE] = Data.points;
  uint(&m_count)[MAX_NODE] = Data.count;
  ulong(&m_dis)[MAX_NODE] = Data.dis;
  double(&m_ans)[MAX_NODE] = Data.ans;
  double(&m_g)[MAX_NODE] = Data.g;
  uint &m_pointNum = Data.pointNum;

  // 找路
  std::priority_queue<Node> pq;
  pq.push(Node{start, 0});
  m_dis[start] = 0;
  m_count[start] = 1;
  while (!pq.empty()) {
    Node head = pq.top();
    pq.pop();
    uint u = head.idx;
    if (head.val > m_dis[u]) continue;
    m_points[++m_pointNum] = u;  //拓扑序
    const DFSEdge *e = &GHead[Head[u]];
    const auto &l = Head[u], &r = Head[u + 1];
    for (uint i = l; i < r; ++i, ++e) {
      const auto &v = e->idx;
      const ulong &w = e->w;
      const ulong &newdis = m_dis[u] + w;
      if (newdis < m_dis[v]) {
        m_count[v] = m_count[u];
        m_dis[v] = newdis;
        pq.push(Node{v, m_dis[v]});
      } else if (newdis == m_dis[v]) {
        m_count[v] += m_count[u];  // s到e.idx的最短路条数
      }
    }
  }

  GetAnswer(Data, start);
  Clear(Data);
}

/*
 * Team: 孤芳自赏
 * No1. chier
 * No2. XDUls
 * No3. yangzhi__
 *
 * [Stragety: Dijkstra with heap手写堆]
 *
 * 这里是一个Dijkstra算法，并用了二叉堆优化，为了看起来好看，我必须加这个注释
 * 1. 手写二叉堆 优先队列记录count
 * 2. 全地球人公用一个更新答案的接口
 * 3. 全地球人公用一个clear的接口
 *
 * 注释不够多，不美观
 */

void Update(SolverData &Data, const uint &x) {
  for (uint i = x, j = x >> 1; j; i = j, j >>= 1) {
    if (Data.dis[Data.head[i]] >= Data.dis[Data.head[j]]) return;
    std::swap(Data.head[i], Data.head[j]);
    std::swap(Data.id[Data.head[i]], Data.id[Data.head[j]]);
  }
}
void Push(SolverData &Data, const uint &x) {
  Data.head[++Data.head[0]] = x;
  Data.id[x] = Data.head[0];
  Update(Data, Data.head[0]);
}
void Pop(SolverData &Data) {
  Data.id[Data.head[1]] = 0;
  Data.id[Data.head[Data.head[0]]] = 1;
  Data.head[1] = Data.head[Data.head[0]--];
  for (uint i = 1, j = 2; j <= Data.head[0]; i = j, j <<= 1) {
    if (Data.dis[Data.head[j + 1]] < Data.dis[Data.head[j]]) ++j;
    if (Data.dis[Data.head[i]] <= Data.dis[Data.head[j]]) return;
    std::swap(Data.head[i], Data.head[j]);
    std::swap(Data.id[Data.head[i]], Data.id[Data.head[j]]);
  }
}

void DijkstraAtHeap(SolverData &Data, const uint &start) {
  uint(&m_points)[MAX_NODE] = Data.points;
  uint(&m_count)[MAX_NODE] = Data.count;
  uint(&m_head)[MAX_NODE] = Data.head;
  uint(&m_id)[MAX_NODE] = Data.id;
  ulong(&m_dis)[MAX_NODE] = Data.dis;
  double(&m_ans)[MAX_NODE] = Data.ans;
  double(&m_g)[MAX_NODE] = Data.g;
  uint &m_pointNum = Data.pointNum;

  // 找路
  Push(Data, start);
  m_dis[start] = 0;
  m_count[start] = 1;
  while (m_head[0]) {
    const auto u = m_head[1];
    Pop(Data);
    m_points[++m_pointNum] = u;  //拓扑序
    const DFSEdge *e = &GHead[Head[u]];
    for (uint i = Head[u]; i < Head[u + 1]; ++i, ++e) {
      const auto &v = e->idx;
      const auto &w = e->w;
      const auto &newdis = m_dis[u] + w;
      if (newdis < m_dis[v]) {
        m_count[v] = m_count[u];
        m_dis[v] = newdis;
        if (!m_id[v]) {
          Push(Data, v);
        } else {
          Update(Data, m_id[v]);
        }
      } else if (newdis == m_dis[v]) {
        m_count[v] += m_count[u];  // s到e.idx的最短路条数
      }
    }
  }

  GetAnswer(Data, start);
  Clear(Data);
}

/*
 * Team: 孤芳自赏
 * No1. chier
 * No2. XDUls
 * No3. yangzhi__
 *
 * [Stragety: SPFA]
 *
 * 这里是一个朴素SPFA算法，为了看起来好看，我必须加这个注释
 * 1. 找最短路
 * 2. 构造拓扑序
 * 3. 全地球人公用一个更新答案的接口
 * 4. 全地球人公用一个clear的接口
 *
 * 注释不够多，不美观
 */

void SPFA(SolverData &Data, const uint &start) {
  uint(&m_points)[MAX_NODE] = Data.points;
  uint(&m_count)[MAX_NODE] = Data.count;
  uint(&m_viscount)[MAX_NODE] = Data.viscount;
  bool(&m_vis)[MAX_NODE] = Data.vis;
  ulong(&m_dis)[MAX_NODE] = Data.dis;
  double(&m_ans)[MAX_NODE] = Data.ans;
  double(&m_g)[MAX_NODE] = Data.g;
  uint &m_pointNum = Data.pointNum;

  std::queue<uint> q;
  q.push(start);
  m_vis[start] = true;
  m_dis[start] = 0;

  while (!q.empty()) {
    const uint u = q.front();
    q.pop();
    m_vis[u] = false;

    const DFSEdge *e = &GHead[Head[u]];
    for (uint i = Head[u]; i < Head[u + 1]; ++i, ++e) {
      const auto &v = e->idx;
      const auto &w = e->w;
      const auto &newdis = m_dis[u] + w;
      if (newdis < m_dis[v]) {
        m_dis[v] = newdis;
        m_viscount[v] = 1;
        if (!m_vis[v]) {
          m_vis[v] = true;
          q.push(v);
        }
      } else if (newdis == m_dis[v]) {
        ++m_viscount[v];
      }
    }
  }

  // 构造拓扑序
  uint l = 0;
  m_points[++m_pointNum] = start;
  m_count[start] = 1;
  while (m_pointNum > l) {
    int u = m_points[++l];
    for (uint i = Head[u]; i < Head[u + 1]; ++i) {
      const auto &e = GHead[i];
      if (m_dis[u] + e.w == m_dis[e.idx]) {
        if (!--m_viscount[e.idx]) m_points[++m_pointNum] = e.idx;
        m_count[e.idx] += m_count[u];
      }
    }
  }

  // 更新答案
  GetAnswer(Data, start);
  Clear(Data);
}

/*
 * Team: 孤芳自赏
 * No1. chier
 * No2. XDUls
 * No3. yangzhi__
 *
 * [Stragety: SPFA + SLF]
 *
 * 这里是一个SPFA + SLF算法，为了看起来好看，我必须加这个注释
 * 1. 常规队列替换为双端队列, 对于一个要加入的点u
 *  - 如果dis[v] < dis[q.front()]，u就放入队首
 *  - 否则放入队尾
 * 2. 构造拓扑序
 * 3. 全地球人公用一个更新答案的接口
 * 4. 全地球人公用一个clear的接口
 *
 * 注释不够多，不美观
 */
void SPFAAtSLF(SolverData &Data, const uint &start) {
  uint(&m_points)[MAX_NODE] = Data.points;
  uint(&m_count)[MAX_NODE] = Data.count;
  uint(&m_viscount)[MAX_NODE] = Data.viscount;
  bool(&m_vis)[MAX_NODE] = Data.vis;
  ulong(&m_dis)[MAX_NODE] = Data.dis;
  double(&m_ans)[MAX_NODE] = Data.ans;
  double(&m_g)[MAX_NODE] = Data.g;
  uint &m_pointNum = Data.pointNum;

  std::deque<uint> q;
  q.push_back(start);
  m_vis[start] = true;
  m_dis[start] = 0;

  while (!q.empty()) {
    const uint u = q.front();
    q.pop_front();
    m_vis[u] = false;

    const DFSEdge *e = &GHead[Head[u]];
    for (uint i = Head[u]; i < Head[u + 1]; ++i, ++e) {
      const auto &v = e->idx;
      const auto &w = e->w;
      const auto &newdis = m_dis[u] + w;
      if (newdis < m_dis[v]) {
        m_dis[v] = newdis;
        m_viscount[v] = 1;
        if (!m_vis[v]) {
          m_vis[v] = true;
          // SLF
          if (q.empty() || m_dis[v] > m_dis[q.front()]) {
            q.push_back(v);
          } else {
            q.push_front(v);
          }
        }
      } else if (newdis == m_dis[v]) {
        ++m_viscount[v];
      }
    }
  }

  // 构造拓扑序
  uint l = 0;
  m_points[++m_pointNum] = start;
  m_count[start] = 1;
  while (m_pointNum > l) {
    int u = m_points[++l];
    for (uint i = Head[u]; i < Head[u + 1]; ++i) {
      const auto &e = GHead[i];
      if (m_dis[u] + e.w == m_dis[e.idx]) {
        if (!--m_viscount[e.idx]) m_points[++m_pointNum] = e.idx;
        m_count[e.idx] += m_count[u];
      }
    }
  }

  // 更新答案
  GetAnswer(Data, start);
  Clear(Data);
}

/*
 * Team: 孤芳自赏
 * No1. chier
 * No2. XDUls
 * No3. yangzhi__
 *
 * [Stragety: SPFA + SLF + LLL]
 *
 * 这里是一个SPFA + SLF + LLL算法，为了看起来好看，我必须加这个注释
 * 1. 对每个要出队的元素u，比较dis[u]和队列中所有dis值的平均值
 *  - 如果dis[u]大，那么将它弹出放到队尾
 *  - 取队首元素再重复判断，直达存在dis[u]小于平均值
 * 2. 构造拓扑序
 * 3. 全地球人公用一个更新答案的接口
 * 4. 全地球人公用一个clear的接口
 *
 * 注释不够多，不美观
 */
void SPFAAtSLFAndLLL(SolverData &Data, const uint &start) {
  uint(&m_points)[MAX_NODE] = Data.points;
  uint(&m_count)[MAX_NODE] = Data.count;
  uint(&m_viscount)[MAX_NODE] = Data.viscount;
  bool(&m_vis)[MAX_NODE] = Data.vis;
  ulong(&m_dis)[MAX_NODE] = Data.dis;
  double(&m_ans)[MAX_NODE] = Data.ans;
  double(&m_g)[MAX_NODE] = Data.g;
  uint &m_pointNum = Data.pointNum;

  std::deque<uint> q;
  q.push_back(start);
  m_vis[start] = true;
  m_dis[start] = 0;

  ulong sum = 0;
  uint cnt = 1;

  while (!q.empty()) {
    uint u = q.front();

    // LLL
    while (cnt * m_dis[u] > sum) {
      q.pop_front();
      q.push_back(u);
      u = q.front();
    }

    q.pop_front();
    --cnt;
    sum -= m_dis[u];
    m_vis[u] = false;

    const DFSEdge *e = &GHead[Head[u]];
    for (uint i = Head[u]; i < Head[u + 1]; ++i, ++e) {
      const auto &v = e->idx;
      const auto &w = e->w;
      const auto &newdis = m_dis[u] + w;
      if (newdis < m_dis[v]) {
        m_dis[v] = newdis;
        m_viscount[v] = 1;
        if (!m_vis[v]) {
          m_vis[v] = true;
          ++cnt;
          sum += newdis;
          // SLF
          if (q.empty() || newdis > m_dis[q.front()]) {
            q.push_back(v);
          } else {
            q.push_front(v);
          }
        }
      } else if (newdis == m_dis[v]) {
        ++m_viscount[v];
      }
    }
  }

  // 构造拓扑序
  uint l = 0;
  m_points[++m_pointNum] = start;
  m_count[start] = 1;
  while (m_pointNum > l) {
    int u = m_points[++l];
    for (uint i = Head[u]; i < Head[u + 1]; ++i) {
      const auto &e = GHead[i];
      if (m_dis[u] + e.w == m_dis[e.idx]) {
        if (!--m_viscount[e.idx]) m_points[++m_pointNum] = e.idx;
        m_count[e.idx] += m_count[u];
      }
    }
  }

  // 更新答案
  GetAnswer(Data, start);
  Clear(Data);
}

/*
 * Team: 孤芳自赏
 * No1. chier
 * No2. XDUls
 * No3. yangzhi__
 *
 * [这里是我Simulation部分]
 *
 * 这里是一个simulation，调度任务，为了看起来好看，我必须加这个注释
 * 1. TopSort: 拓扑排序，为了gay掉可以删除的点
 * 2. GetJob: atomic维护任务队列，动态调度任务
 * 3. FindTask: 启动线程，选择哪种策略开始Gay
 *
 * 注释不够多，不美观
 */

void TopSort() {
  std::queue<uint> q;
  for (uint i = 0; i < g_NodeNum; ++i) {
    if (BackLen[i] == 0) q.push(i);
  }
  while (!q.empty()) {
    uint u = q.front();
    q.pop();

    Top[u] = true;

    for (uint i = Head[u]; i < Head[u + 1]; ++i) {
      const auto &e = GHead[i];
      if (HeadLen[u] == 1) {
        Label[e.idx] += (Label[u] + 1);
      }
      if (--BackLen[e.idx] <= 0) q.push(e.idx);
    }
  }
}

inline void GetJob(uint &job) {
  static uint l = 0;
  static std::atomic_flag lock = ATOMIC_FLAG_INIT;
  while (lock.test_and_set())
    ;
  job = l < g_NodeNum ? l++ : -1;
  lock.clear();
}

void printProcess(const uint &job) {
  if (job % 1000 == 0) {
    std::cerr << "[" << job << "/" << g_NodeNum << "]\n";
  }
}

void FindTask(uint pid) {
  auto &Data = SovData[pid];
  for (uint i = 0; i < g_NodeNum; ++i) {
    Data.count[i] = 0;
    Data.ans[i] = 0;
    Data.dis[i] = UINT64_MAX;
    Data.g[i] = 0;
    Data.head[i] = 0;
    Data.id[i] = 0;
    Data.vis[i] = false;
    Data.viscount[i] = 0;
  }

  uint job = 0;
  while (true) {
    GetJob(job);
    if (job == -1) break;
    if (Top[job] && HeadLen[job] == 1) continue;

    /*
     * GAY B, GAY wo
     *
     * IfSparseGraph = true 表示稀疏图
     *
     * 1. Dijkstra() 稀疏图
     * 2. DijkstraAtHeap() 稠密图
     * 3. SPFA() 我也不知道能咋
     * 4. SPFAAtSLF() magic优化
     * 5. SPFAAtSLFAndLLL() magic 优化
     */

    // Dijkstra(Data, job);
    // DijkstraAtHeap(Data, job);
    // SPFA(Data, job);
    // SPFAAtSLF(Data, job);
    SPFAAtSLFAndLLL(Data, job);

#ifdef DEBUG
    printProcess(job);
#endif
  }
}

/*
 * Team: 孤芳自赏
 * No1. chier
 * No2. XDUls
 * No3. yangzhi__
 *
 * [这里是我处理答案部分]
 *
 * 这里是一个SaveAnswer，为了看起来好看，我必须加这个注释
 * 1. Find: 启动T个线程,完了之后合并答案
 * 2. SaveAnswer: 答案排一下顺序,然后保存到文件
 *
 * 注释不够多，不美观
 */

typedef std::pair<uint, double> prud;
std::vector<prud> Answer;  // 最终结果

void Find() {
  std::thread Th[T];
  for (uint i = 0; i < T; ++i) Th[i] = std::thread(FindTask, i);
  for (auto &it : Th) it.join();

  Answer = std::vector<prud>(g_NodeNum);
  for (uint i = 0; i < g_NodeNum; ++i) {
    Answer[i].first = IDDom[i];
    Answer[i].second = 0;
    for (const auto &data : SovData) {
      Answer[i].second += data.ans[i];
    }
  }
}

/************************** Algorithm End ****************************/

void SaveAnswer() {
  std::sort(Answer.begin(), Answer.end(), [](const prud &p1, const prud &p2) {
    if (p1.second == p2.second) return p1.first < p2.first;
    return p1.second > p2.second;
  });

  uint up = g_NodeNum > 100 ? 100 : g_NodeNum;

  char *buffer = new char[100 * 100];
  char *ptr = buffer;

  for (uint i = 0; i < up; ++i) {
    const auto &it = Answer[i];
    ptr += sprintf(ptr, "%d,%.3lf\n", it.first, it.second);
  }

  FILE *fp = fopen(RESULT, "w");
  fwrite(buffer, 1, ptr - buffer, fp);
  fclose(fp);
}

/*
 * Team: 孤芳自赏
 * No1. chier
 * No2. XDUls
 * No3. yangzhi__
 *
 * [这里是我分析地图特征的部分]
 *
 * 这里是一个分析地图特征，为了看起来好看，我必须加这个注释
 * 1. UnionSet: 并查集瞅一瞅有多少个联通快
 * 2. AnalysisGraph: 分析一下稠密图还是稀疏图，打印一下地图情况
 *
 * 注释不够多，不美观
 */

void UnionSet() {
  for (uint i = 0; i <= g_NodeNum; ++i) Fa[i] = i;
  for (uint i = 0; i < g_EdgeNum; ++i) {
    auto &e = Edges[i];
    uint x = getFa(e.u), y = getFa(e.v);
    Fa[x] = Fa[y];
  }

  for (uint i = 0; i < g_NodeNum; ++i) {
    uint x = getFa(i);
    if (x == i) BlockNum++;
    ++Blocks[x].first;
    Blocks[x].second += HeadLen[i];
  }
}

void AnalysisGraph() {
#ifdef DEBUG
  Timer t;
#endif
  // E <= V * 10 -> 稀疏图
  IfSparseGraph = g_EdgeNum <= 10 * g_NodeNum ? true : false;
  UnionSet();
  uint cnt = 0;
  for (uint i = 0; i < g_NodeNum; ++i) {
    if (Top[i] && HeadLen[i] == 1) {
      ++cnt;
    }
  }
#ifdef DEBUG
  Color::green();
  std::cerr << "==================================\n";
  std::cerr << "* 地图: " << (IfSparseGraph ? "稀疏图" : "稠密图") << "\n";
  std::cerr << "* 结点: " << g_NodeNum << "\n";
  std::cerr << "* 边数: " << g_EdgeNum << "\n";
  std::cerr << "* 集合: " << BlockNum << "\n";
  std::cerr << "* 删除: " << cnt << "\n";
  std::cerr << "* cost: " << t.elapsed() << "s\n";
  std::cerr << "==================================\n";
  Color::reset();
#endif
}

/*
 * Team: 孤芳自赏
 * No1. chier
 * No2. XDUls
 * No3. yangzhi__
 *
 * [这里是程序main函数]，这里要是看不懂就别写代码了
 *
 * 注释不够多，不美观
 */
int main() {
  std::cerr << std::fixed << std::setprecision(4);

  LoadData();
  TopSort();
  AnalysisGraph();
  Find();
  SaveAnswer();

  return 0;
}