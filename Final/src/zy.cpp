#include <bits/stdc++.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstring>
#include <ext/pb_ds/priority_queue.hpp>
#include <iostream>
#include <map>
#include <queue>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;
#ifdef LOCAL

#include <sys/time.h>

//#define TESTFILE "/data/onlinedata/data6/test_data.txt"
//#define TESTFILE "/data/1w/test_data.txt"
//#define TESTFILE "data/newsmall/node1wedge2w.txt"
#define TESTFILE "../data/std1/test_data.txt"
//#define TESTFILE "data/oneedge.txt"
#define RESULT "../data/std1/result.txt"
#else
#define TESTFILE "/data/test_data.txt"
#define RESULT "/projects/student/result.txt"
#endif

#define uint uint32_t
#define ulong uint64_t

/**全局常量**/
const uint kMaxE = 3000024;    //最大边数
const uint kMaxN = 5000024;    //最大点数
const uint kThreadNum = 8;     //进程/线程数量
const uint kMaxIdStrLen = 12;  //整数表示的ID有多长
enum { SPFA, DIJ, DIJHEAP, DIJSEGTREE };
const int GLOBAL_STRATEGY = DIJSEGTREE;
/**结构体**/
// typedef pair<ulong, uint> Pli;
// typedef __gnu_pbds::priority_queue<Pli, greater<Pli>,
// __gnu_pbds::pairing_heap_tag> Heap;

class HashMap {
  //优化的unordered map from github
 private:
  const static uint max_size = 1900037;

 public:
  uint count = 0;
  uint hash_vals[max_size];
  struct Bucket {
    uint key;
    int val = -1;
  } buckets[max_size];

 private:
  uint hash(const uint &i, const int &k) {
    return (k % max_size + i * (4516783 - k % 4516783)) % max_size;
  }

 public:
  int operator[](const uint &key) {
    for (int i = 0; i < max_size; ++i) {
      const uint hash_val = hash(key, i);
      if (buckets[hash_val].val == -1)
        return -1;
      else if (buckets[hash_val].key == key)
        return buckets[hash_val].val;
    }
    return -1;
  }

  void insert(const uint &key, const uint &value) {
    for (int i = 0; i < max_size; ++i) {
      const uint hashval = hash(key, i);
      if (buckets[hashval].val == -1) {
        buckets[hashval].key = key;
        buckets[hashval].val = value;
        hash_vals[count++] = hashval;
        return;
      } else if (buckets[hashval].key == key)
        return;
    }
  }

  bool tryInsert(const uint &key, const uint &value) {
    for (int i = 0; i < max_size; ++i) {
      const uint hashval = hash(key, i);
      if (buckets[hashval].val == -1) {
        buckets[hashval].key = key;
        buckets[hashval].val = value;
        hash_vals[count++] = hashval;
        return true;
      } else if (buckets[hashval].key == key)
        return false;
    }
    return false;
  }
};

template <class T>
class Vec {
  uint elem_cnt;
  uint elem_capacity;
  T *elem;

 public:
  /** Default constructor */
  explicit Vec<T>(uint s)
      : elem_cnt(0),
        elem_capacity(s),
        elem{static_cast<T *>(::operator new(sizeof(T) * elem_capacity))} {};
  explicit Vec<T>() : elem_cnt(0), elem_capacity(0), elem(nullptr){};
  Vec<T>(const Vec<T> &other)
      : elem_cnt(other.elem_cnt),
        elem_capacity(other.elem_capacity),
        elem(new T[elem_cnt]) {
    for (uint i = 0; i < elem_cnt; ++i) elem[i] = other.elem[i];
  }

  Vec<T> &operator=(const Vec<T> &other) {
    for (uint i = 0; i < elem_cnt; ++i) elem[i] = other.elem[i];
    elem_cnt = other.elem_cnt;
    elem_capacity = other.elem_capacity;
    return *this;
  }

  T *begin() { return elem; }

  T *begin() const { return elem; }

  T *end() { return (elem + elem_cnt); }

  T *end() const { return (elem + elem_cnt); }

  T &operator[](uint n) { return elem[n]; }

  T &operator[](uint n) const { return elem[n]; }

  inline uint size() const { return elem_cnt; }

  template <typename... Args>
  inline void emplace_back(Args &&... args) noexcept {
    if (elem_cnt == elem_capacity) {
      reserve(elem_capacity * 2 + 1);
    }
    new (&elem[elem_cnt++]) T(std::forward<Args>(args)...);
  }
  void reserve(const uint &n) {
    /*        if (elem_capacity < n) {
                T* tmp=static_cast<T*>(realloc(elem,n* sizeof(T)));
                elem=tmp;
                elem_capacity = n;
            }*/
    if (elem_capacity < n) {
      elem_capacity = n;
      T *tmp{static_cast<T *>(::operator new(sizeof(T) * elem_capacity))};
      memmove(tmp, &elem[0], sizeof(T) * elem_cnt);
      ::operator delete(elem);
      elem = tmp;
    }
  }
  inline bool empty() { return elem_cnt == 0; }
  inline void push_back(T &value) {
    if (elem_cnt >= elem_capacity) {
      reserve(elem_capacity * 2 + 1);
    }
    elem[elem_cnt++] = value;
  }
  inline void push_back(const T &value) {
    if (elem_cnt >= elem_capacity) {
      reserve(elem_capacity * 2 + 1);
    }
    elem[elem_cnt++] = value;
  }
  inline void clear() { elem_cnt = 0; }
  //只有在第二次改元素才能用
  inline void reInit(const T &value) {
    if (elem_cnt >= elem_capacity) {
      reserve(elem_capacity * 2 + 1);
    }
    elem_cnt = 1;
    elem[0] = value;
  }
};

struct RawEdge {
  uint u;
  uint v;
  ulong w;
};
struct BriefEdge {
  uint v;
  ulong w;
};
struct IdString {
  uint len;
  char val[kMaxIdStrLen];
};
struct SPPQNode {
  uint u;
  ulong dist;

  bool operator<(const SPPQNode &x) const { return dist > x.dist; }
};
struct SegmentTree {
  uint n = 0;
  uint offset = 1;
  vector<ulong> min_w;
  // vector<uint> min_pos;
  SegmentTree() = default;
  explicit SegmentTree(uint node_num) {
    n = node_num;
    min_w = vector<ulong>(4 * node_num + 4, UINT64_MAX);
    while (offset < n) offset = offset << 1;
  }

  void modify(uint pos, ulong w) {
    // update(1, 0, offset - 1, pos, w);
    update2(pos, w);
  }

  void backUpdate(uint at) {
    min_w[at] = min_w[at << 1] < min_w[at << 1 | 1] ? min_w[at << 1]
                                                    : min_w[at << 1 | 1];
    if (at > 1) backUpdate(at >> 1);
  }
  void update2(uint pos, ulong w) {
    min_w[pos + offset] = w;
    backUpdate((pos + offset) >> 1);
  }

  uint query(uint at, uint left, uint right) {
    if (left == right) return left;
    uint mid = (left + right) >> 1;
    if (min_w[at << 1] < min_w[at << 1 | 1])
      return query(at << 1, left, mid);
    else
      return query(at << 1 | 1, mid + 1, right);
  }

  uint top() { return query(1, 0, offset - 1); }
};

struct WorkerLoadInfo {
  uint node_num = 0;
  uint edge_num = 0;
  uint map_offset = 0;
  // unordered_map<uint, uint> id_map;
  HashMap id_map;
  RawEdge raw_edge[kMaxE];
  uint mod_edge_num[kThreadNum]{};
  RawEdge mod_edge[kThreadNum][kMaxE / kThreadNum + 1];
  uint mod_y_num[kThreadNum]{};
  uint mod_y[kThreadNum][kMaxE / kThreadNum + 1];
};
// TODO 连通分量  稠密时用指针++访问
struct WorkerFindInfo {
  ulong dist[kMaxN];

  Vec<uint> prev[kMaxN];
  Vec<uint> back[kMaxN];
  double res[kMaxN];
  bool vis[kMaxN];
  uint path_num[kMaxN];
  uint tp[kMaxN];  //拓扑

  SegmentTree seg_tree;
  vector<double> result;
};

/**全局变量*/
uint g_node_num = 0;
uint g_edge_num = 0;

uint g_out_list[kMaxN];
BriefEdge g_edge[kMaxE];
IdString g_node_str[kMaxN];
uint g_unorder_node_id[kMaxN];
RawEdge g_raw_edge[kMaxE];
WorkerLoadInfo w_load_info[kThreadNum];
WorkerFindInfo w_find_info[kThreadNum];

vector<double> g_result;  //关键中心性

atomic_flag g_find_lock = ATOMIC_FLAG_INIT;
uint g_find_num = 0;

uint g_in_list_size[kMaxN];
uint g_in_list_offset[kMaxN];
uint g_in_list[kMaxN];
BriefEdge g_reverse_edge[kMaxE];
/**读数据*/
void addEdge(uint x, uint y, ulong w, WorkerLoadInfo &data) {
  if (w == 0) return;
  uint mod_x_pid = x % kThreadNum;
  RawEdge &e = data.mod_edge[mod_x_pid][data.mod_edge_num[mod_x_pid]];
  e.u = x;
  e.v = y;
  e.w = w;
  data.mod_edge_num[mod_x_pid]++;
  uint mod_y_pid = y % kThreadNum;
  data.mod_y[mod_y_pid][data.mod_y_num[mod_y_pid]] = y;
  data.mod_y_num[mod_y_pid]++;
}

void analyzeBuffer(const char *buffer, uint buffer_size, uint pid) {
  const char *p = buffer;
  uint x, y;
  ulong w;
  while (p < buffer + buffer_size) {
    x = 0;
    y = 0;
    w = 0;
    while (*p != ',') {
      x = (x << 3) + (x << 1) + *p - '0';
      ++p;
    }
    ++p;
    while (*p != ',') {
      y = (y << 3) + (y << 1) + *p - '0';
      ++p;
    }
    ++p;
    while (*p != '\r' && *p != '\n') {
      w = (w << 3) + (w << 1) + *p - '0';
      ++p;
    }
    if (*p == '\r') ++p;
    ++p;
    addEdge(x, y, w, w_load_info[pid]);
  }
}

void readBuffer() {
  uint fd = open(TESTFILE, O_RDONLY);
  uint buffer_size = lseek(fd, 0, SEEK_END);
  char *file_buffer =
      (char *)mmap(NULL, buffer_size, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);

  thread pool[kThreadNum - 1];
  uint step = buffer_size / kThreadNum;
  char *p = file_buffer;
  char *q = file_buffer + step;
  for (uint i = 0; i < kThreadNum - 1; ++i) {
    if (q - file_buffer >= buffer_size) {
      pool[i] = thread(analyzeBuffer, p, 0, i);
      continue;
    }
    while ((*q) != '\n') {
      q++;
    }
    pool[i] = thread(analyzeBuffer, p, q - p + 1, i);
    p = ++q;
    q = p + step;
  }
  analyzeBuffer(p, file_buffer + buffer_size - p, kThreadNum - 1);
  for (auto &th : pool) th.join();
}

void sortRawEdge(uint pid) {
  WorkerLoadInfo &data = w_load_info[pid];
  RawEdge *p = data.raw_edge;
  data.edge_num = 0;
  for (uint i = 0; i < kThreadNum; ++i) {
    memcpy(p, w_load_info[i].mod_edge[pid],
           w_load_info[i].mod_edge_num[pid] * sizeof(RawEdge));
    p += w_load_info[i].mod_edge_num[pid];
    data.edge_num += w_load_info[i].mod_edge_num[pid];
  }
  // sort(data.raw_edge, data.raw_edge + data.edge_num, cmpEdge);
  sort(data.raw_edge, data.raw_edge + data.edge_num,
       [&](const RawEdge &a, const RawEdge &b) {
         return a.u != b.u ? (a.u < b.u) : (a.v < b.v);
       });
}

void establishIdMap(uint pid) {
  WorkerLoadInfo &data = w_load_info[pid];
  HashMap &id_map = data.id_map;
  uint edge_num = data.edge_num;
  uint node_num = 0;
  uint last = -1;
  for (uint i = 0; i < edge_num; ++i) {
    RawEdge &e = data.raw_edge[i];
    if (e.u != last) {
      last = e.u;
      id_map.insert(last, node_num++);
    }
  }

  for (uint i = 0; i < kThreadNum; ++i) {
    for (uint j = 0; j < w_load_info[i].mod_y_num[pid]; ++j) {
      if (id_map.tryInsert(w_load_info[i].mod_y[pid][j], node_num)) node_num++;
    }
  }

  data.node_num = node_num;
}

void calculateMapOffset() {
  w_load_info[0].map_offset = 0;
  g_node_num = w_load_info[0].node_num;
  for (uint i = 0; i < kThreadNum - 1; ++i) {
    w_load_info[i + 1].map_offset =
        w_load_info[i].map_offset + w_load_info[i].node_num;
    g_node_num += w_load_info[i + 1].node_num;
  }
}

void mergeUniqueId(uint pid) {
  const uint offset = w_load_info[pid].map_offset;
  HashMap &id_map = w_load_info[pid].id_map;
  for (uint i = 0; i < id_map.count; ++i) {
    auto &p = id_map.buckets[id_map.hash_vals[i]];
    p.val += offset;
    g_unorder_node_id[p.val] = p.key;
  }
}

uint *idArgSort() {
  uint *arg_list = new uint[g_node_num];
  uint *arg_reflect = new uint[g_node_num];
  for (uint i = 0; i < g_node_num; ++i) arg_list[i] = i;

  std::sort(arg_list, arg_list + g_node_num, [](int pos1, int pos2) {
    return (g_unorder_node_id[pos1] < g_unorder_node_id[pos2]);
  });
  for (uint i = 0; i < g_node_num; ++i) {
    arg_reflect[arg_list[i]] = i;
  }
  return arg_reflect;
}

void intToStr(uint x, IdString &str) {
  if (x == 0) {
    str.val[0] = '0';
    str.val[1] = ',';
    str.len = 2;
  } else {
    char tmp[kMaxIdStrLen];
    uint idx = kMaxIdStrLen;
    tmp[--idx] = ',';
    while (x > 0) {
      tmp[--idx] = x % 10 + '0';
      x /= 10;
    }
    memcpy(str.val, tmp + idx, kMaxIdStrLen - idx);
    str.len = kMaxIdStrLen - idx;
  }
}

void mergeIdmap(const uint *arg_reflect, uint pid) {
  HashMap &id_map = w_load_info[pid].id_map;
  for (uint i = 0; i < id_map.count; ++i) {
    auto &p = id_map.buckets[id_map.hash_vals[i]];
    p.val = arg_reflect[p.val];
    intToStr(p.key, g_node_str[p.val]);
  }
}

void mapEdgeUV(uint pid) {
  WorkerLoadInfo &data = w_load_info[pid];
  for (uint i = 0; i < data.edge_num; ++i) {
    RawEdge &e = data.raw_edge[i];
    // uint x_id = data.id_map[e.u];
    uint x_id = data.id_map[e.u];
    e.u = x_id;
    uint y_pid = e.v % kThreadNum;

    auto y_id = w_load_info[y_pid].id_map[e.v];
    if (y_id == -1) {
      printf("cannot map %d\n", e.v);
    }
    e.v = y_id;
  }
}

struct PQNode {
  uint first;
  uint second;

  bool operator<(const PQNode &b)
      const  //写在里面只用一个b，但是要用const和&修饰，并且外面还要const修饰;
  {
    return first > b.first;
  }
};

//多线程比单线程慢
void mergeRawEdge() {
  g_edge_num = 0;
  priority_queue<PQNode> edge_heap;
  uint head[kThreadNum] = {0};
  for (uint i = 0; i < kThreadNum; ++i) {
    if (w_load_info[i].edge_num > 0) {
      edge_heap.push({w_load_info[i].raw_edge[0].u, i});
      // cout<<w_load_info[i].raw_edge[0].u<<"
      // "<<w_load_info[i].raw_edge[0].v<<" "<<w_load_info[i].raw_edge[0].w<<"in
      // "<<i<<endl;
    }
  }
  uint next;
  uint last_u;
  RawEdge *e_master = g_raw_edge;
  while (!edge_heap.empty()) {
    next = edge_heap.top().second;
    edge_heap.pop();

    RawEdge *e_worker = &w_load_info[next].raw_edge[head[next]];
    last_u = e_worker->u;
    // printf("edge %d %d %d\n",e_worker->u,e_worker->v,e_worker->w);
    uint len = 0;
    while (last_u == e_worker->u &&
           e_worker - w_load_info[next].raw_edge < w_load_info[next].edge_num) {
      ++e_worker;
      ++len;
      // printf("edge %d %d %d\n",e_worker->u,e_worker->v,e_worker->w);
    }
    // printf("len %d lastu %d\n",len,last_u);
    memcpy(e_master, &w_load_info[next].raw_edge[head[next]],
           len * sizeof(RawEdge));
    e_master += len;

    head[next] += len;
    if (head[next] < w_load_info[next].edge_num)
      edge_heap.push({w_load_info[next].raw_edge[head[next]].u, next});
  }
  g_edge_num = e_master - g_raw_edge;
  // printf("edge num %d\n",g_edge_num);
}

void buildGraph(uint pid) {
  uint u_last = 0;
  RawEdge *e = g_raw_edge;
  for (uint i = 0; i < g_edge_num; ++i) {
    if (e->u % kThreadNum == pid) {
      if (e->u != u_last) {
        g_out_list[e->u] = i;
      }
      u_last = e->u;
    }
    if (e->v % kThreadNum == pid) {
      // g_in_list_old[e->v].emplace_back(e->u, e->w);
      g_in_list_size[e->v]++;
      BriefEdge &b = g_edge[i];
      b.v = e->v;
      b.w = e->w;
    }
    ++e;
  }
  g_out_list[g_node_num] = g_edge_num;  //必须加 否则最后一项不对
}

void buildGraphFixTail() {
  g_out_list[g_node_num] = g_edge_num;
  uint tail = g_node_num - 1;
  uint first_unzero = g_raw_edge[0].u;

  for (uint i = tail; i > first_unzero; --i) {
    if (g_out_list[i] == 0) {
      g_out_list[i] = g_out_list[i + 1];
    }
  }
}
void accumulateInListSize() {
  for (uint i = 1; i < g_node_num; ++i) {  // g_inlist[0]=0
    g_in_list[i] = g_in_list[i - 1] + g_in_list_size[i - 1];
    g_in_list_offset[i] = g_in_list[i];
  }
  g_in_list[g_node_num] = g_edge_num;
}

void setInList(uint pid) {
  for (int i = g_edge_num - 1; i >= 0; --i) {
    RawEdge &e = g_raw_edge[i];
    if (e.v % kThreadNum == pid) {
      g_reverse_edge[g_in_list_offset[e.v]].v = e.u;
      g_reverse_edge[g_in_list_offset[e.v]].w = e.w;
      ++g_in_list_offset[e.v];
    }
  }
}

void loadData() {
  // stage0
  readBuffer();

  // stage1
  thread pool[kThreadNum - 1];
  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(sortRawEdge, i + 1);
  }
  sortRawEdge(0);
  for (auto &th : pool) th.join();
  // stage2

  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(establishIdMap, i + 1);
  }
  establishIdMap(0);
  for (auto &th : pool) th.join();

  // stage3
  calculateMapOffset();
  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(mergeUniqueId, i + 1);
  }
  mergeUniqueId(0);
  for (auto &th : pool) th.join();

  // stage4
  uint *arg_reflect = idArgSort();
  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(mergeIdmap, arg_reflect, i + 1);
  }
  mergeIdmap(arg_reflect, 0);
  for (auto &th : pool) th.join();

  // stage5
  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(mapEdgeUV, i + 1);
  }
  mapEdgeUV(0);
  for (auto &th : pool) th.join();

  // stage6
  mergeRawEdge();

  for (uint i = 0; i < g_edge_num; ++i) {
    // printf("%d %d %d\n",g_raw_edge[i].u,g_raw_edge[i].v,g_raw_edge[i].w);
  }
  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(buildGraph, i + 1);
  }
  buildGraph(0);
  for (auto &th : pool) th.join();

  buildGraphFixTail();
  accumulateInListSize();
  for (uint i = 0; i < kThreadNum - 1; ++i) pool[i] = thread(setInList, i);
  setInList(kThreadNum - 1);
  for (auto &th : pool) th.join();

  printf("N:%d E:%d\n", g_node_num, g_edge_num);
}

// TODO 转账金额为0怎么处理
// TODO SPFA优化
// TODO dijkstra
/**1.SPFA naive*/

/**2.dij naive*/

// TODO 循环展开
uint g_one_in[kThreadNum]{0};
void findAtDijkstra(uint s, WorkerFindInfo &info, uint pid) {
  ulong(&dist)[kMaxN] = info.dist;
  Vec<uint>(&prev)[kMaxN] = info.prev;
  Vec<uint>(&back)[kMaxN] = info.back;
  double(&res)[kMaxN] = info.res;
  bool(&vis)[kMaxN] = info.vis;

  uint(&path_num)[kMaxN] = info.path_num;
  uint(&tp)[kMaxN] = info.tp;
  int rank = 0;
  priority_queue<SPPQNode> q;
  dist[s] = 0;
  q.push({s, 0});
  path_num[s] = 1;

  while (!q.empty()) {
    const uint u = q.top().u;
    q.pop();
    if (vis[u]) continue;
    vis[u] = true;
    tp[rank++] = u;
    for (uint i = g_out_list[u]; i < g_out_list[u + 1]; ++i) {
      const BriefEdge &e = g_edge[i];
      const uint v = e.v;
      // const uint w = e.w;
      ulong new_dist = dist[u] + e.w;
      if (new_dist < dist[v]) {
        dist[v] = new_dist;
        prev[v].reInit(u);
        q.push({v, new_dist});
      } else if (new_dist == dist[v]) {
        prev[v].push_back(u);
      }
    }
  }

  for (uint i = 1; i < rank; ++i) {
    // if (i == s) continue;
    for (const auto &v : prev[tp[i]]) {
      back[v].push_back(tp[i]);
      path_num[tp[i]] += path_num[v];
    }
  }

  for (int i = rank - 1; i >= 1; i--) {
    for (const auto &v : back[tp[i]]) {
      res[tp[i]] +=
          (res[v] + double(1)) * double(path_num[tp[i]]) / double(path_num[v]);
    }
  }
  for (uint i = 1; i < rank; ++i) {
    info.result[tp[i]] += res[tp[i]];
  }

  for (uint i = g_in_list[s]; i < g_in_list[s + 1]; ++i) {
    uint v = g_reverse_edge[i].v;
    if (g_out_list[v + 1] - g_out_list[v] == 1 && !vis[v]) {
      g_one_in[pid]++;
    }
  }

  for (uint r = 0; r < rank; ++r) {
    uint i = tp[r];
    dist[i] = UINT64_MAX;
    prev[i].clear();
    back[i].clear();
    res[i] = 0;
    vis[i] = false;
    path_num[i] = 0;
  }
}

/**3.dij 配对堆*/
// TODO 配对堆  斐波那契堆  拓扑排序

// TODO memset
void findAtDijSegmentTree(uint s, WorkerFindInfo &info) {
  ulong(&dist)[kMaxN] = info.dist;
  Vec<uint>(&prev)[kMaxN] = info.prev;
  Vec<uint>(&back)[kMaxN] = info.back;
  double(&res)[kMaxN] = info.res;
  uint(&path_num)[kMaxN] = info.path_num;
  uint(&tp)[kMaxN] = info.tp;
  int rank = 0;

  dist[s] = 0;
  path_num[s] = 1;
  SegmentTree &q = info.seg_tree;
  q.modify(s, 0);

  while (q.min_w[1] != UINT64_MAX) {
    const uint u = q.top();
    q.modify(u, UINT64_MAX);
    tp[rank++] = u;
    for (uint i = g_out_list[u]; i < g_out_list[u + 1]; ++i) {
      const BriefEdge &e = g_edge[i];
      const uint v = e.v;
      const ulong w = e.w;
      ulong new_dist = dist[u] + w;
      if (new_dist < dist[v]) {
        dist[v] = new_dist;
        prev[v].reInit(u);
        q.modify(v, new_dist);
      } else if (new_dist == dist[v]) {
        prev[v].push_back(u);
      }
    }
  }

  for (uint i = 1; i < rank; ++i) {
    // if (i == s) continue;
    for (const auto &v : prev[tp[i]]) {
      back[v].push_back(tp[i]);
      path_num[tp[i]] += path_num[v];
    }
  }

  for (int i = rank - 1; i >= 1; i--) {
    for (const auto &v : back[tp[i]]) {
      res[tp[i]] +=
          (res[v] + double(1)) * double(path_num[tp[i]]) / double(path_num[v]);
    }
  }
  for (uint i = 1; i < rank; ++i) {
    info.result[tp[i]] += res[tp[i]];
  }

  for (uint r = 0; r < rank; ++r) {
    uint i = tp[r];
    dist[i] = UINT64_MAX;
    prev[i].clear();
    back[i].clear();
    res[i] = 0;
    path_num[i] = 0;
  }
}

void initInfo(WorkerFindInfo &info) {
  ulong(&dist)[kMaxN] = info.dist;
  Vec<uint>(&prev)[kMaxN] = info.prev;
  Vec<uint>(&back)[kMaxN] = info.back;
  double(&res)[kMaxN] = info.res;
  bool(&vis)[kMaxN] = info.vis;

  uint(&path_num)[kMaxN] = info.path_num;
  uint(&tp)[kMaxN] = info.tp;
  int rank = 0;

  for (uint i = 0; i < g_node_num; ++i) {
    dist[i] = UINT64_MAX;
  }
  for (uint i = 0; i < g_node_num; ++i) {
    prev[i].clear();
    back[i].clear();
  }
  for (uint i = 0; i < g_node_num; ++i) {
    res[i] = 0;
  }
  for (uint i = 0; i < g_node_num; ++i) {
    vis[i] = false;
  }
  for (uint i = 0; i < g_node_num; ++i) {
    path_num[i] = 0;
  }
}

void findShortestPath(uint pid) {
  w_find_info[pid].result = vector<double>(g_node_num, 0);
  initInfo(w_find_info[pid]);
  if (GLOBAL_STRATEGY == DIJSEGTREE) {
    w_find_info[pid].seg_tree = SegmentTree(g_node_num);
  }

  uint next_node;
  while (true) {
    while (g_find_lock.test_and_set()) {
    }
    next_node = g_find_num < g_node_num ? g_find_num++ : g_node_num;
    g_find_lock.clear();
    if (next_node >= g_node_num) break;
#ifdef LOCAL
    if (next_node % 1000 == 0) printf("%d\n", next_node);
#endif
    switch (GLOBAL_STRATEGY) {
      // case SPFA:findAtSPFA(next_node, w_find_info[pid]);break;
      case DIJ:
        findAtDijkstra(next_node, w_find_info[pid], pid);
        break;
        // case DIJHEAP:findAtDijPairHeap(next_node, w_find_info[pid]);break;
      case DIJSEGTREE:
        findAtDijSegmentTree(next_node, w_find_info[pid]);
        break;
      default:
        break;
    }
  }
}

void findPath() {
  thread pool[kThreadNum - 1];
  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(findShortestPath, i);
  }
  findShortestPath(kThreadNum - 1);
  for (auto &th : pool) th.join();
}

void writeResult() {
  g_result = vector<double>(g_node_num);
  for (uint i = 0; i < kThreadNum; ++i) {
    for (uint j = 0; j < g_node_num; ++j) {
      g_result[j] += w_find_info[i].result[j];
    }
  }

  // uint *arg_list = new uint[g_node_num];
  vector<uint> arg_list(g_node_num);
  for (uint i = 0; i < g_node_num; ++i) arg_list[i] = i;

  std::sort(arg_list.begin(), arg_list.end(), [&](int pos1, int pos2) {
    return abs(g_result[pos1] - g_result[pos2]) > 0.0001
               ? g_result[pos1] > g_result[pos2]
               : pos1 < pos2;
  });

  FILE *fd = fopen(RESULT, "w");
  char *answer_buffer = new char[50 * 100];
  char *p = answer_buffer;
  uint answer_size = g_node_num < 100 ? g_node_num : 100;

  for (uint i = 0; i < answer_size; ++i) {
    memcpy(p, g_node_str[arg_list[i]].val,
           g_node_str[arg_list[i]].len);  // TODO 前面不用存了 这里直接int转就行
    p += g_node_str[arg_list[i]].len;
    p += sprintf(p, "%.3lf\n", g_result[arg_list[i]]);
  }
  fwrite(answer_buffer, p - answer_buffer, 1, fd);
}

void solve() {
#ifdef LOCAL
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif

  loadData();

#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("load data time %fs\n", t2 - t1);
#endif
  findPath();
#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t3 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("find cycle time %fs\n", t3 - t2);
#endif
  writeResult();
}

int main() { solve(); }