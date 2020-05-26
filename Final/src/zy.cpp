#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <map>
#include <queue>
#include <stack>
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
#define TESTFILE "data/practice/test_data.txt"
#define RESULT "/tmp/result.txt"
#else
#define TESTFILE "/data/test_data.txt"
#define RESULT "/projects/student/result.txt"
#endif

#define uint uint32_t
#define ulong uint64_t

/**全局常量**/
const uint kMaxE = 3000024;    //最大边数
const uint kMaxN = 3000024;    //最大点数
const uint kThreadNum = 8;     //进程/线程数量
const uint kMaxIdStrLen = 12;  //整数表示的ID有多长
enum { SPFA, DIJSPARSE, DIJSEGTREE };
int GLOBAL_STRATEGY = DIJSPARSE;

/**结构体**/
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
  uint idx;
  uint father;
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
    uint x = pos + offset;
    min_w[x] = w;
    // backUpdate((pos + offset) >> 1);
    while (x > 1) {
      x >>= 1;
      min_w[x] =
          min_w[x << 1] < min_w[x << 1 | 1] ? min_w[x << 1] : min_w[x << 1 | 1];
    }
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
  char vis[kMaxN];
  uint path_num[kMaxN];
  uint tp[kMaxN];  //拓扑
  uint n_cut[kMaxN];
  ulong dist[kMaxN];
  double res[kMaxN];
  double result[kMaxN];
  priority_queue<SPPQNode> pq;
  vector<uint> in_all;
};

/**全局变量*/
uint g_node_num = 0;
uint g_edge_num = 0;

uint g_out_list[kMaxN];
BriefEdge g_edge[kMaxE];
uint g_in_list[kMaxN];
BriefEdge g_reverse_edge[kMaxE];

IdString g_node_str[kMaxN];
uint g_unorder_node_id[kMaxN];
RawEdge g_raw_edge[kMaxE];
uint g_in_list_size[kMaxN];
uint g_in_list_offset[kMaxN];
WorkerLoadInfo w_load_info[kThreadNum];
WorkerFindInfo w_find_info[kThreadNum];

vector<double> g_result;  //关键中心性
atomic_flag g_find_lock = ATOMIC_FLAG_INIT;
uint g_find_num = 0;

bool g_visit[kMaxN];
vector<uint> g_arg_insize_order;

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
  sort(data.raw_edge, data.raw_edge + data.edge_num,
       [&](const RawEdge &a, const RawEdge &b) {
         return a.u != b.u ? (a.u < b.u) : (a.w < b.w);
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

//多线程比单线程慢
void mergeRawEdge() {
  struct PQNode {
    uint first;
    uint second;

    bool operator<(const PQNode &b)
        const  //写在里面只用一个b，但是要用const和&修饰，并且外面还要const修饰;
    {
      return first > b.first;
    }
  };
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

void getInSizeOrder() {
  g_arg_insize_order = vector<uint>(g_node_num);
  for (uint i = 0; i < g_node_num; i++) g_arg_insize_order[i] = i;
  sort(g_arg_insize_order.begin(), g_arg_insize_order.end(),
       [&](uint pos1, uint pos2) {
         return g_in_list[pos1 + 1] - g_in_list[pos1] >
                g_in_list[pos2 + 1] - g_in_list[pos2];
       });
  // g_in_list_size[pos1]>g_in_list[pos2];
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

  getInSizeOrder();
  printf("N:%d E:%d\n", g_node_num, g_edge_num);
}

/**根据稀疏稠密选方法*/
// TODO 提前算outlistsize=1
void dijNoCut(uint s, WorkerFindInfo &info) {
  char(&vis)[kMaxN] = info.vis;
  ulong(&dist)[kMaxN] = info.dist;
  uint(&path_num)[kMaxN] = info.path_num;
  uint(&tp)[kMaxN] = info.tp;
  double(&res)[kMaxN] = info.res;
  double(&result)[kMaxN] = info.result;
  int rank = 0;
  priority_queue<SPPQNode> &q = info.pq;
  // priority_queue<SPPQNode> q;
  dist[s] = 0;
  q.push({s, s, 0});
  path_num[s] = 1;
  stack<pair<uint, uint>> stk;
  while (!q.empty()) {
    const uint u = q.top().idx;
    const uint father = q.top().father;
    const ulong dist_u = q.top().dist;
    q.pop();
    if (!vis[u]) {
      vis[u] = 1;
      tp[rank++] = u;
      stk.push({father, u});
      const BriefEdge *e = &g_edge[g_out_list[u]];
      for (uint i = g_out_list[u]; i < g_out_list[u + 1]; ++i, ++e) {
        const uint v = e->v;
        const ulong w = e->w;
        if (!vis[v]) {
          ulong new_dist = dist_u + w;
          if (new_dist < dist[v]) {
            dist[v] = new_dist;
            path_num[v] = path_num[u];
            q.push({v, u, new_dist});
          } else if (new_dist == dist[v]) {
            path_num[v] += path_num[u];
            q.push({v, u, new_dist});
          }
        }
      }
    } else if (dist_u == dist[u]) {
      stk.push({father, u});
    }
  }
  while (stk.size() > 1) {
    const uint u = stk.top().first;
    const uint v = stk.top().second;
    stk.pop();
    res[u] += (res[v] + double(1)) * double(path_num[u]) / double(path_num[v]);
  }
  for (uint i = 1; i < rank; ++i) {
    result[tp[i]] += res[tp[i]];
  }
  for (uint r = 0; r < rank; ++r) {
    uint i = tp[r];
    dist[i] = UINT64_MAX;
    res[i] = 0;
    vis[i] = false;
    path_num[i] = 0;
  }
}

uint handle(uint s, WorkerFindInfo &info) {
  vector<uint> search_list;
  for (uint i = g_in_list[s]; i < g_in_list[s + 1]; ++i) {
    const uint v = g_reverse_edge[i].v;
    if (g_out_list[v + 1] - g_out_list[v] == 1 && !info.vis[v]) {
      while (g_find_lock.test_and_set()) {
      }
      if (!g_visit[v]) {
        g_visit[v] = true;
        g_find_lock.clear();
        search_list.push_back(v);
      } else
        g_find_lock.clear();
    }
  }
  double new_res = info.res[s] + 1.0;
  for (const auto v : search_list) {
    info.res[v] = new_res;
    info.n_cut[s] += handle(v, info);
  }

  info.in_all.push_back(s);
  return info.n_cut[s] + 1;
}

void calResultDijSparse(uint s, uint rank, WorkerFindInfo &info) {
  ulong(&dist)[kMaxN] = info.dist;
  double(&res)[kMaxN] = info.res;
  char(&vis)[kMaxN] = info.vis;
  uint(&path_num)[kMaxN] = info.path_num;
  uint(&tp)[kMaxN] = info.tp;
  uint(&n_cut)[kMaxN] = info.n_cut;
  vector<uint> &in_all = info.in_all;
  handle(s, info);

  for (uint i = 1; i < rank; ++i) {
    info.result[tp[i]] += res[tp[i]] * (1 + n_cut[s]);
  }
  for (const auto &v : in_all) {
    info.result[v] += res[v] * n_cut[v];
  }
  for (uint r = 0; r < rank; ++r) {
    uint i = tp[r];
    dist[i] = UINT64_MAX;
    res[i] = 0;
    vis[i] = 0;
    path_num[i] = 0;
  }
  for (const auto v : in_all) {
    n_cut[v] = 0;
    res[v] = 0;
  }
  in_all.clear();
}

// TODO 循环展开
void dijkstraSparse(uint s, WorkerFindInfo &info) {
  char(&vis)[kMaxN] = info.vis;
  ulong(&dist)[kMaxN] = info.dist;
  uint(&path_num)[kMaxN] = info.path_num;
  uint(&tp)[kMaxN] = info.tp;
  double(&res)[kMaxN] = info.res;
  double(&result)[kMaxN] = info.result;
  int rank = 0;
  priority_queue<SPPQNode> &q = info.pq;
  // priority_queue<SPPQNode> q;
  dist[s] = 0;
  q.push({s, s, 0});
  path_num[s] = 1;
  stack<pair<uint, uint>> stk;
  while (!q.empty()) {
    const uint u = q.top().idx;
    const uint father = q.top().father;
    const ulong dist_u = q.top().dist;
    q.pop();
    if (!vis[u]) {
      vis[u] = 1;
      tp[rank++] = u;
      stk.push({father, u});
      const BriefEdge *e = &g_edge[g_out_list[u]];
      for (uint i = g_out_list[u]; i < g_out_list[u + 1]; ++i, ++e) {
        const uint v = e->v;
        const ulong w = e->w;
        if (!vis[v]) {
          ulong new_dist = dist_u + w;
          if (new_dist < dist[v]) {
            dist[v] = new_dist;
            path_num[v] = path_num[u];
            q.push({v, u, new_dist});
          } else if (new_dist == dist[v]) {
            path_num[v] += path_num[u];
            q.push({v, u, new_dist});
          }
        }
      }
    } else if (dist_u == dist[u]) {
      stk.push({father, u});
    }
  }
  while (stk.size() > 1) {
    const uint u = stk.top().first;
    const uint v = stk.top().second;
    stk.pop();
    res[u] += (res[v] + double(1)) * double(path_num[u]) / double(path_num[v]);
  }
  calResultDijSparse(s, rank, info);
  //    for (uint i = 1; i < rank; ++i) {
  //        result[tp[i]] += res[tp[i]];
  //    }
  //    for (uint r = 0; r < rank; ++r) {
  //        uint i = tp[r];
  //        dist[i] = UINT64_MAX;
  //        res[i] = 0;
  //        vis[i] = false;
  //        path_num[i] = 0;
  //    }
}

void initInfo(WorkerFindInfo &info) {
  ulong(&dist)[kMaxN] = info.dist;
  double(&res)[kMaxN] = info.res;
  char(&vis)[kMaxN] = info.vis;
  uint(&path_num)[kMaxN] = info.path_num;
  uint(&tp)[kMaxN] = info.tp;
  double(&result)[kMaxN] = info.result;
  uint(&n_cut)[kMaxN] = info.n_cut;
  for (uint i = 0; i < g_node_num; ++i) {
    dist[i] = UINT64_MAX;
    res[i] = 0;
    vis[i] = false;
    path_num[i] = 0;
    result[i] = 0;
    result[i] = 0;
    n_cut[i] = 0;
  }
}

void findShortestPath(uint pid) {
  initInfo(w_find_info[pid]);
  uint next_node;
  while (true) {
    while (g_find_lock.test_and_set()) {
    }
    next_node =
        g_find_num < g_node_num ? g_arg_insize_order[g_find_num++] : g_node_num;
    if (g_visit[next_node]) {
      g_find_lock.clear();
      if (next_node >= g_node_num) break;
      continue;
    }
    g_visit[next_node] = true;
    g_find_lock.clear();
    if (next_node >= g_node_num) break;
#ifdef LOCAL
    if (g_find_num % 1000 == 0) printf("%d\n", g_find_num);
#endif
    switch (GLOBAL_STRATEGY) {
      case DIJSPARSE:
        dijkstraSparse(next_node, w_find_info[pid]);
        break;
      case SPFA:
        dijNoCut(next_node, w_find_info[pid]);
        break;
      default:
        break;
    }
    /*        if (g_find_num > 6000)
                exit(0);*/
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

int partition(vector<uint> &a, int l, int r) {
  swap(a[l], a[l + rand() % (r - l + 1)]);  //随机化
  int j = l;
  double pivot = g_result[a[l]];
  int pos = a[l];
  for (int i = l + 1; i <= r; i++) {
    if (g_result[a[i]] - pivot > 0.0001 ||
        (abs(g_result[a[i]] - pivot) <= 0.0001 && a[i] < pos))
      swap(a[i], a[++j]);
  }
  swap(a[l], a[j]);
  return j;
}

void select(vector<uint> &a, int l, int r, int k) {
  int j = partition(a, l, r);
  int k1 = k - 1;
  if (k1 == j)
    return;
  else if (k1 < j)
    select(a, l, j - 1, k);
  else
    select(a, j + 1, r, k);
}

vector<uint> getLeastNumbers(vector<uint> &arr, int k) {
  if (k != 0) select(arr, 0, arr.size() - 1, k);
  arr.resize(k);
  return arr;
}

void writeResult() {
  g_result = vector<double>(g_node_num);
  for (uint i = 0; i < kThreadNum; ++i) {
    for (uint j = 0; j < g_node_num; ++j) {
      g_result[j] += w_find_info[i].result[j];
    }
  }

  uint answer_size = g_node_num < 100 ? g_node_num : 100;

  vector<uint> arg_list(g_node_num);
  for (uint i = 0; i < g_node_num; ++i) arg_list[i] = i;

  getLeastNumbers(arg_list, answer_size);

  std::sort(arg_list.begin(), arg_list.end(), [&](int pos1, int pos2) {
    return abs(g_result[pos1] - g_result[pos2]) > 0.0001
               ? g_result[pos1] > g_result[pos2]
               : pos1 < pos2;
  });

  FILE *fd = fopen(RESULT, "w");
  char *answer_buffer = new char[50 * 100];
  char *p = answer_buffer;

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

  // setStrategy();

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