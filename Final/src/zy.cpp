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
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;
#ifdef LOCAL

#include <sys/time.h>
#define TESTFILE "/data/tiny/test_data.txt"
//#define TESTFILE "/data/1w/test_data.txt"
#define RESULT "/tmp/result.txt"
#else
#define TESTFILE "/data/test_data.txt"
#define RESULT "/projects/student/result.txt"
#endif

#define uint uint32_t
#define ulong uint64_t

/**全局常量**/
const uint kMaxE = 3000024;    //最大边数
const uint kMaxN = 1000024;    //最大点数
const uint kThreadNum = 8;     //进程/线程数量
const uint kMaxIdStrLen = 12;  //整数表示的ID有多长

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

struct WorkerLoadInfo {
  uint node_num = 0;
  uint edge_num = 0;
  uint map_offset = 0;
  // unordered_map<uint, uint> id_map;
  HashMap id_map;
  RawEdge raw_edge[kMaxE];
  uint mod_edge_num[kThreadNum];
  RawEdge mod_edge[kThreadNum][kMaxE / kThreadNum + 1];
  uint mod_y_num[kThreadNum];
  uint mod_y[kThreadNum][kMaxN / kThreadNum + 1];
};
struct WorkerFindInfo {
  vector<bool> found;
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

/**读数据*/
void addEdge(uint x, uint y, ulong w, WorkerLoadInfo &data) {
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
  cout << "x1\n";

  for (uint i = 0; i < g_edge_num; ++i) {
    // printf("%d %d %d\n",g_raw_edge[i].u,g_raw_edge[i].v,g_raw_edge[i].w);
  }
  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(buildGraph, i + 1);
  }
  buildGraph(0);
  for (auto &th : pool) th.join();

  buildGraphFixTail();

#ifdef LOCAL
  printf("N:%d E:%d\n", g_node_num, g_edge_num);
#endif
}

// TODO 转账金额为0怎么处理
// TODO SPFA优化
// TODO dijkstra
void findAtSPFA(uint s, WorkerFindInfo &info) {
  vector<ulong> dist = vector<ulong>(g_node_num, UINT64_MAX);
  dist[s] = 0;
  vector<vector<uint>> prev = vector<vector<uint>>(g_node_num);
  vector<vector<uint>> back = vector<vector<uint>>(g_node_num);
  vector<double> res = vector<double>(g_node_num, 0);
  queue<uint> q;  //搜索队列
  q.push(s);
  vector<bool> in_queue = vector<bool>(g_node_num, false);
  in_queue[s] = true;
  vector<uint> path_num = vector<uint>(g_node_num, 0);
  path_num[s] = 1;

  while (!q.empty()) {
    const uint u = q.front();
    q.pop();
    in_queue[u] = false;
    for (uint i = g_out_list[u]; i < g_out_list[u + 1]; ++i) {
      const auto &e = g_edge[i];
      if (dist[u] + e.w > dist[e.v]) continue;
      if (dist[u] + e.w == dist[e.v]) {
        prev[e.v].push_back(u);
      } else {
        dist[e.v] = dist[u] + e.w;
        prev[e.v] = {u};
        if (!in_queue[e.v]) {
          q.push(e.v);
          in_queue[e.v] = true;
        }
      }
    }
  }

  for (uint i = 0; i < g_node_num; ++i) {
    if (i == s) continue;
    for (const auto &v : prev[i]) {
      back[v].push_back(i);
    }
  }
  vector<uint> vis_num(g_node_num, 0);
  queue<uint> tmp;

  tmp.push(s);
  while (!tmp.empty()) {
    const uint u = tmp.front();
    tmp.pop();
    if (vis_num[u] < prev[u].size()) {
      tmp.push(u);
      continue;
    }
    for (const auto &v : back[u]) {
      path_num[v] += path_num[u];
      if (vis_num[v] == 0) tmp.push(v);
      vis_num[v]++;
    }
  }

  for (uint i = 0; i < g_node_num; ++i) {
    vis_num[i] = 0;
  }

  for (uint i = 0; i < g_node_num; ++i) {
    if (i == s) continue;
    if (back[i].size() == 0) tmp.push(i);
  }
  while (!tmp.empty()) {
    const uint u = tmp.front();
    tmp.pop();
    if (vis_num[u] < back[u].size()) {
      tmp.push(u);
      continue;
    }
    for (const auto &v : prev[u]) {
      if (v == s) continue;
      res[v] +=
          (res[u] + double(1)) * double(path_num[v]) / double(path_num[u]);
      if (vis_num[v] == 0) tmp.push(v);
      vis_num[v]++;
    }
  }
  for (uint i = 0; i < g_node_num; ++i) {
    info.result[i] += res[i];
  }

  // TODO 分两段可以少一个判断
}

struct PQNode2 {
  uint u;
  ulong dist;

  bool operator<(const PQNode2 &x) const { return dist > x.dist; }
};

// TODO 用spfa的inqueue优化dij可行吗
void findAtDijkstra(uint s, WorkerFindInfo &info) {
  vector<ulong> dist = vector<ulong>(g_node_num, UINT64_MAX);
  dist[s] = 0;
  vector<vector<uint>> prev = vector<vector<uint>>(g_node_num);
  vector<vector<uint>> back = vector<vector<uint>>(g_node_num);
  vector<double> res = vector<double>(g_node_num, 0);
  vector<bool> vis(g_node_num, false);
  priority_queue<PQNode2> q;
  q.push({s, 0});
  vector<uint> path_num = vector<uint>(g_node_num, 0);
  path_num[s] = 1;

  while (!q.empty()) {
    const uint u = q.top().u;
    q.pop();
    if (vis[u]) continue;
    vis[u] = true;
    for (uint i = g_out_list[u]; i < g_out_list[u + 1]; ++i) {
      const BriefEdge &e = g_edge[i];
      const uint v = e.v;
      const ulong w = e.w;
      ulong new_dist = dist[u] + w;
      if (new_dist < dist[v]) {
        dist[v] = new_dist;
        prev[v] = {u};
        q.push({v, new_dist});
      } else if (new_dist == dist[v]) {
        prev[v].push_back(u);
      }
    }
  }
  for (uint i = 0; i < g_node_num; ++i) {
    if (i == s) continue;
    for (const auto &v : prev[i]) {
      back[v].push_back(i);
    }
  }
  vector<uint> vis_num(g_node_num, 0);
  queue<uint> tmp;

  tmp.push(s);
  while (!tmp.empty()) {
    const uint u = tmp.front();
    tmp.pop();
    if (vis_num[u] < prev[u].size()) {
      tmp.push(u);
      continue;
    }
    for (const auto &v : back[u]) {
      path_num[v] += path_num[u];
      if (vis_num[v] == 0) tmp.push(v);
      vis_num[v]++;
    }
  }

  for (uint i = 0; i < g_node_num; ++i) {
    vis_num[i] = 0;
  }

  for (uint i = 0; i < g_node_num; ++i) {
    if (i == s) continue;
    if (back[i].size() == 0) tmp.push(i);
  }
  while (!tmp.empty()) {
    const uint u = tmp.front();
    tmp.pop();
    if (vis_num[u] < back[u].size()) {
      tmp.push(u);
      continue;
    }
    for (const auto &v : prev[u]) {
      if (v == s) continue;
      res[v] +=
          (res[u] + double(1)) * double(path_num[v]) / double(path_num[u]);
      if (vis_num[v] == 0) tmp.push(v);
      vis_num[v]++;
    }
  }
  for (uint i = 0; i < g_node_num; ++i) {
    info.result[i] += res[i];
  }
}

void findShortestPath(uint pid) {
  w_find_info[pid].result = vector<double>(g_node_num, 0);
  uint next_node;
  while (true) {
    while (g_find_lock.test_and_set()) {
    }
    next_node = g_find_num < g_node_num ? g_find_num++ : g_node_num;
    g_find_lock.clear();
    if (next_node >= g_node_num) break;
    // if (next_node % 1000 == 0)
    // printf("%d\n", next_node);
    // findAtSPFA(next_node, w_find_info[pid]);
    findAtDijkstra(next_node, w_find_info[pid]);
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
    return g_result[pos1] != g_result[pos2] ? g_result[pos1] > g_result[pos2]
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