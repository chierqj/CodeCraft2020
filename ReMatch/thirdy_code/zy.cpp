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
#define TESTFILE "../data/18875018/test_data.txt"
//#define TESTFILE "testdata/sfn.txt"
//#define TESTFILE "testdata/test_data_massive.txt"
#define RESULT "/dev/shm/result.txt"
#else
#define TESTFILE "/data/test_data.txt"
#define RESULT "/projects/student/result.txt"
#endif

#define uint uint32_t

/**全局常量**/
const uint kMaxResult = 22000128;  //最大结果数
const uint kMaxE = 3000024;        //最大边数
const uint kMaxN = 2000024;        //最大点数
const uint kThreadNum = 4;         //进程/线程数量
const uint kModVal = kThreadNum - 1;
const uint kMaxIdStrLen = 12;         //整数表示的ID有多长
const uint kUselessFlg = 2147123456;  //标记没有出边的点
const uint kBufferNum = 64;           //
const uint kMULT3MAX = 715827882;
const uint kMULT5MAX = 429496729;

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
};

struct RawEdge {
  uint u;
  uint v;
  uint w;
};
// typedef pair<uint,uint> uintE;
typedef uint uintE;
struct IdString {
  uint len;
  char val[kMaxIdStrLen];
};
struct Result {
  vector<uint> path[5];
};
struct WorkerLoadInfo {
  uint node_num = 0;
  uint edge_num = 0;
  uint map_offset = 0;
  // unordered_map<uint, uint> id_map;
  uint mod_edge_num[kThreadNum];
  RawEdge mod_edge[kThreadNum][kMaxE / kThreadNum];
  HashMap id_map;
  RawEdge raw_edge[kMaxE];
};
struct WorkerFindInfo {
  uint visit_num = 0;      //记录反向搜索改了多少点的dist
  uint visit_list[kMaxN];  //改了哪些点的dist
  char dist[kMaxN];        //反向搜索得到的点离st的距离
  uint end_weight
      [kMaxN];  //记录反向搜索第一层的点，用于比较最后一条边和第一条边的weight是否合法
};
struct WorkerWriteInfo {
  uint start_idx = 0;
  uint end_idx = 0;
  uint answer_len = 0;
  uint is_finish = 0;
  char *answer_buffer;  //[kMaxResult * 7 * kMaxIdStrLen / kBufferNum];
};

/**全局变量*/
// graph&find
uint g_node_num = 0;
uint g_edge_num = 0;

uint g_in_list_size[kMaxN];
uint g_in_list_offset[kMaxN];
uint g_in_list[kMaxN];
uint g_out_list[kMaxN];
uintE g_edge[kMaxE][2];
uintE g_reverse_edge[kMaxE][2];
vector<Result> g_result;
// read temp
IdString g_node_str[kMaxN];
uint g_unorder_node_id[kMaxN];
RawEdge g_raw_edge[kMaxE];
// V4.3 将worker变量全部打包
// worker
WorkerLoadInfo w_load_info[kThreadNum];
WorkerFindInfo w_find_info[kThreadNum];
WorkerWriteInfo w_write_info[kBufferNum];
// lock
atomic_flag g_find_lock = ATOMIC_FLAG_INIT;
uint g_find_num = 0;
atomic_flag g_buffer_lock = ATOMIC_FLAG_INIT;
uint g_prepared_buffer_num = 0;

uint result_num = 0;
// atomic_flag g_write_lock = ATOMIC_FLAG_INIT;
uint g_next_buffer = 0;
vector<uint> g_flower_point;

/**读数据*/
/** 一个新的逻辑，如果一个点只出现在y中而不出现在x中，那它必然不构成任何环，因此我们只需要为rawedge中的u建立映射即可
 * 又由于边是排序过的，因此仅当遍历中发现不一样的u时才尝试插入到idmap
 *
 * TODO 两种str策略，一种是加在rawedge里，另一种是只存intid，后面再反解析成str
 * */
void addEdge(uint x, uint y, uint w, WorkerLoadInfo &info) {
  uint mod_pid = x % kThreadNum;
  RawEdge &e = info.mod_edge[mod_pid][info.mod_edge_num[mod_pid]];
  e.u = x;
  e.v = y;
  e.w = w;
  info.mod_edge_num[mod_pid]++;
}

void analyzeBuffer(const char *buffer, uint buffer_size, uint pid) {
  const char *p = buffer;
  uint x, y, w;
  /*    while (q - buffer < max_buffer_size) {
          while (*q != split)
              q++;

          parseInteger(p, q, x);
          p = ++q;
          while (*q != split)
              q++;

          parseInteger(p, q, y);
          p = ++q;

          while (*q != '\r' && *q != '\n')
              q++;
          parseInteger(p, q, w);
          if (*q == '\r')
              ++q;
          p = ++q;
          addEdge(x, y, w, pid);
      }*/
  while (p < buffer + buffer_size) {
    x = 0;
    y = 0;
    w = 0;
    while (*p != ',') {
      x = (*p - '0') + (x << 3) + (x << 1);
      ++p;
    }
    ++p;
    while (*p != ',') {
      y = (*p - '0') + (y << 3) + (y << 1);
      ++p;
    }
    ++p;
    while (*p != '\r' && *p != '\n') {
      w = (*p - '0') + (w << 3) + (w << 1);
      ++p;
    }
    if (*p == '\r') ++p;
    ++p;
    addEdge(x, y, w, w_load_info[pid]);
  }
}

void readBuffer() {
  /*    FILE *fd = fopen(TESTFILE, "r");
      char *file_buffer;
      fseek(fd, 0, SEEK_END);//将文件内部的指针指向文件末尾
      uint buffer_size =
     ftell(fd);//获取文件长度，（得到文件位置指针当前位置相对于文件首的偏移字节数）
      rewind(fd);//将文件内部的指针重新指向一个流的开头
      file_buffer = (char *) malloc(buffer_size * sizeof(char) +
     1);//申请内存空间，lsize*sizeof(char)是为了更严谨，16位上char占一个字符，其他机器上可能变化
      fread(file_buffer, 1, buffer_size, fd);//将pfile中内容读入pread指向内存中
      file_buffer[buffer_size] = '\0';*/
  // mmap 稳定快0.01
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
    q = p + buffer_size;
  }
  analyzeBuffer(p, file_buffer + buffer_size - p, kThreadNum - 1);
  for (auto &th : pool) th.join();
}

void sortRawEdge(uint pid) {
  WorkerLoadInfo &info = w_load_info[pid];
  RawEdge *p = info.raw_edge;
  info.edge_num = 0;
  for (uint i = 0; i < kThreadNum; ++i) {
    memcpy(p, w_load_info[i].mod_edge[pid],
           w_load_info[i].mod_edge_num[pid] * sizeof(RawEdge));
    p += w_load_info[i].mod_edge_num[pid];
    info.edge_num += w_load_info[i].mod_edge_num[pid];
  }
  // sort(info.raw_edge, info.raw_edge + info.edge_num, cmpEdge);
  sort(info.raw_edge, info.raw_edge + info.edge_num,
       [&](const RawEdge &a, const RawEdge &b) {
         return a.u != b.u ? (a.u < b.u) : (a.v < b.v);
       });
}

void establishIdMap(uint pid) {
  WorkerLoadInfo &info = w_load_info[pid];
  // unordered_map<uint, uint> &id_map = info.id_map;
  // id_map.reserve(kMaxN / 4);
  // id_map.rehash(kMaxN / 4);
  HashMap &id_map = info.id_map;
  uint edge_num = info.edge_num;
  uint node_num = 0;
  uint last = -1;
  for (uint i = 0; i < edge_num; ++i) {
    RawEdge &e = info.raw_edge[i];
    if (e.u != last) {
      last = e.u;
      // id_map.insert(make_pair(last, node_num++));
      id_map.insert(last, node_num++);
    }
  }
  info.node_num = node_num;
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
  //    for (auto &p:w_load_info[pid].id_map) {
  //        p.second += offset;
  //        g_unorder_node_id[p.second] = p.first;
  //    }
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

  sort(arg_list, arg_list + g_node_num, [](int pos1, int pos2) {
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
  //    for (auto &p:w_load_info[pid].id_map) {
  //        p.second = arg_reflect[p.second];
  //        intToStr(p.first, g_node_str[p.second]);
  //    }
  HashMap &id_map = w_load_info[pid].id_map;
  for (uint i = 0; i < id_map.count; ++i) {
    auto &p = id_map.buckets[id_map.hash_vals[i]];
    p.val = arg_reflect[p.val];
    intToStr(p.key, g_node_str[p.val]);
  }
}

void mapEdgeUV(uint pid) {
  /*    unordered_map<uint,uint>::iterator map_end[kThreadNum];
      for (uint i=0;i<kThreadNum;++i){
          map_end[i]=g_worker_id_map[i].end();
      }*/
  WorkerLoadInfo &info = w_load_info[pid];
  for (uint i = 0; i < info.edge_num; ++i) {
    RawEdge &e = info.raw_edge[i];
    // uint x_id = info.id_map[e.u];
    uint x_id = info.id_map[e.u];
    e.u = x_id;
    uint y_pid = e.v % kThreadNum;
    uint y_id = kUselessFlg;
    //        auto it = w_load_info[y_pid].id_map.find(e.v);//TODO 提前存end
    //        会变慢 if (it != w_load_info[y_pid].id_map.end()) {
    //            y_id = it->second;
    //        }
    auto it = w_load_info[y_pid].id_map[e.v];
    if (it != -1) {
      y_id = it;
    }
    e.v = y_id;
  }
}

uint getNextRawEdge(uint (&head)[kThreadNum]) {
  uint min_val = 2147483647;
  uint min_idx = kThreadNum;
  for (uint i = 0; i < kThreadNum; ++i) {
    if (w_load_info[i].raw_edge[head[i]].u < min_val &&
        head[i] < w_load_info[i].edge_num) {
      min_val = w_load_info[i].raw_edge[head[i]].u;
      min_idx = i;
    }
  }
  return min_idx;
}

/*
void mergeRawEdge() {
    g_edge_num = 0;
    uint head[kThreadNum] = {0};
    uint next = getNextRawEdge(head);
    uint last_u;
    RawEdge *e_master = g_raw_edge;
    while (next < kThreadNum) {//getNext在都找完后会返回kThreadNum
        RawEdge *e_worker = &w_load_info[next].raw_edge[head[next]];
        last_u = e_worker->u;
        while (last_u == e_worker->u) {
            if (e_worker->v != kUselessFlg) {
                memcpy(e_master, e_worker, sizeof(RawEdge));
                ++e_master;
            }
            ++e_worker;
        }
        head[next] = e_worker - w_load_info[next].raw_edge;
        next = getNextRawEdge(head);
    }
    g_edge_num = e_master - g_raw_edge;
}
*/

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
void mergeRawEdgePQ() {
  g_edge_num = 0;
  priority_queue<PQNode> edge_heap;
  uint head[kThreadNum] = {0};
  for (uint i = 0; i < kThreadNum; ++i) {
    if (w_load_info[i].edge_num > 0) {
      edge_heap.push({w_load_info[i].raw_edge[0].u, i});
      // cout<<w_load_info[i].raw_edge[0].u<<"in  "<<i<<endl;
    }
  }
  uint next;
  uint last_u;
  RawEdge *e_master = g_raw_edge;
  while (!edge_heap.empty()) {
    // cout<<edge_heap.top().first<<" "<<edge_heap.top().second<<endl;
    next = edge_heap.top().second;
    edge_heap.pop();

    RawEdge *e_worker = &w_load_info[next].raw_edge[head[next]];
    last_u = e_worker->u;
    while (last_u == e_worker->u) {
      if (e_worker->v != kUselessFlg) {
        memcpy(e_master, e_worker, sizeof(RawEdge));
        ++e_master;
      }
      ++e_worker;
    }
    head[next] = e_worker - w_load_info[next].raw_edge;
    if (head[next] < w_load_info[next].edge_num)
      edge_heap.push({w_load_info[next].raw_edge[head[next]].u, next});
  }
  g_edge_num = e_master - g_raw_edge;
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
      uintE *b = g_edge[i];
      b[0] = e->v;
      b[1] = e->w;
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
  for (uint i = 1; i < g_node_num + 1; ++i) {  // g_inlist[0]=0
    g_in_list[i] = g_in_list[i - 1] + g_in_list_size[i - 1];
    g_in_list_offset[i] = g_in_list[i];
  }
  // g_in_list[g_node_num] = g_edge_num;
}

void setInList(uint pid) {
  RawEdge *e = &g_raw_edge[0];
  for (uint i = 0; i < g_edge_num; ++i) {
    uint right_v = e->v + 1;
    if (e->v % kThreadNum == pid) {
      g_reverse_edge[g_in_list_offset[right_v] - 1][0] = e->u;
      g_reverse_edge[g_in_list_offset[right_v] - 1][1] = e->w;
      --g_in_list_offset[right_v];
    }
    ++e;
  }
}

/**
 * 读数据的思路：
 * 1.扫一遍所有数据，将原始边保存为互不冲突的T(线程数)组，(按edge.u%4分配)
 * Edge=(u,v,w)
 * 2.各个线程对各自分配到的边进行排序并建立id_map，为了避免冲突建立id_map时只需要对edge.u进行映射即可。
 * 每个线程映射过后，没有被
 *  任何线程映射的点就是没有出边的点，可以忽略
 * 3.各个线程中id_map里的点是unique的，将它们合并，然后argsort得到一个映射，其能保证映射后点的index和id的顺序完全相同。
 *  同时保存它们的str以便后续写答案时memcpy
 * 4.将各个线程中有效的边重新映射，然后合并；保存其简化的副本即只有edge.v和edge.w，以优化找环的访存
 * 5.多线程建图，按前向星存(由于边已经被sort过了，这时前向星中可以只存边的开始和结束序号)
 * 6.结束序号就是下一个点的开始序号，因此多存一个额外假点的开始序号就行
 * 7.反向不再需要sort，由于正向已经完全sort，预先计算好偏移量反向遍历依次插入即可
 * **/
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
  mergeRawEdgePQ();
  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(buildGraph, i + 1);
  }
  buildGraph(0);
  for (auto &th : pool) th.join();
  buildGraphFixTail();

  // stage7
  accumulateInListSize();
  for (uint i = 0; i < kThreadNum - 1; ++i) pool[i] = thread(setInList, i);
  setInList(kThreadNum - 1);
  for (auto &th : pool) th.join();
}

/**找环*/

inline bool cmpWeight(const uint &x, const uint &y) {
  return (y > kMULT5MAX || x <= (y + (y << 2))) &&
         (x > kMULT3MAX || (x + (x << 1)) >= y);
}

void backward3(const uint st, WorkerFindInfo &info) {
  for (int i = 0; i < info.visit_num; ++i) {
    info.dist[info.visit_list[i]] = 0;
  }
  info.visit_list[0] = st;
  info.visit_num = 1;
  info.dist[st] = 7;
  const uintE *e1 = g_reverse_edge[g_in_list[st]];
  for (uint i1 = g_in_list[st]; i1 < g_in_list[st + 1]; ++i1) {
    const uint v1 = *e1;
    ++e1;
    if (v1 <= st) break;
    const uint w1 = *e1;
    ++e1;
    info.dist[v1] = 7;  // 111
    info.end_weight[v1] = w1;
    info.visit_list[info.visit_num++] = v1;
    const uintE *e2 = g_reverse_edge[g_in_list[v1]];
    for (uint i2 = g_in_list[v1]; i2 < g_in_list[v1 + 1]; ++i2) {
      const uint v2 = *e2;
      ++e2;
      if (v2 <= st) break;
      const uint w2 = *e2;
      ++e2;
      if (!cmpWeight(w2, w1)) continue;
      info.dist[v2] |= 3;  // 011
      info.visit_list[info.visit_num++] = v2;
      const uintE *e3 = g_reverse_edge[g_in_list[v2]];
      for (uint i3 = g_in_list[v2]; i3 < g_in_list[v2 + 1]; ++i3) {
        const uint v3 = *e3;
        ++e3;
        if (v3 <= st) break;
        const uint w3 = *e3;  // HINT:位置不能错
        ++e3;
        if (v3 == v1) continue;
        if (!cmpWeight(w3, w2)) continue;
        info.dist[v3] |= 1;  // 001
        info.visit_list[info.visit_num++] = v3;
      }
    }
  }
}

void forward7(const uint st, WorkerFindInfo &info) {
  vector<uint>(&st_result)[5] = g_result[st].path;
  // st_result[4].reserve(400);

  uint i1 = g_out_list[st];
  for (const uintE *e1 = g_edge[i1]; i1 < g_out_list[st + 1]; ++i1) {
    const uint v1 = *e1;
    ++e1;
    const uint w1 = *e1;
    ++e1;
    if (v1 < st) continue;

    uint i2 = g_out_list[v1];
    for (const uintE *e2 = g_edge[i2]; i2 < g_out_list[v1 + 1]; ++i2) {
      const uint v2 = *e2;
      ++e2;
      const uint w2 = *e2;
      ++e2;
      if (v2 <= st || !cmpWeight(w1, w2)) continue;

      uint i3 = g_out_list[v2];
      for (const uintE *e3 = g_edge[i3]; i3 < g_out_list[v2 + 1]; ++i3) {
        const uint v3 = *e3;
        ++e3;
        const uint w3 = *e3;
        ++e3;
        if (v3 < st || v3 == v1 || !cmpWeight(w2, w3)) continue;
        if (v3 == st) {
          if (cmpWeight(w3, w1)) {
            // st_result[0].insert(st_result[0].end(), {st, v1, v2});
            // testnum++;
            // st_result[0].push_back(st);
            st_result[0].push_back(v1);
            st_result[0].push_back(v2);
          }
          continue;
        };
        uint i4 = g_out_list[v3];
        for (const uintE *e4 = g_edge[i4]; i4 < g_out_list[v3 + 1]; ++i4) {
          const uint v4 = *e4;
          ++e4;
          const uint w4 = *e4;
          ++e4;
          if (!(info.dist[v4] & 1) || v4 == v2 || v4 == v1 ||
              !cmpWeight(w3, w4))
            continue;
          else if (v4 == st) {
            if (cmpWeight(w4, w1)) {
              st_result[1].push_back(v1);
              st_result[1].push_back(v2);
              st_result[1].push_back(v3);
            }
            continue;
          }

          uint i5 = g_out_list[v4];
          for (const uintE *e5 = g_edge[i5]; i5 < g_out_list[v4 + 1]; ++i5) {
            const uint v5 = *e5;
            ++e5;
            const uint w5 = *e5;
            ++e5;
            if (!(info.dist[v5] & 2) || v5 == v1 || v5 == v2 || v5 == v3 ||
                !cmpWeight(w4, w5))
              continue;
            else if (v5 == st) {
              if (cmpWeight(w5, w1)) {
                st_result[2].push_back(v1);
                st_result[2].push_back(v2);
                st_result[2].push_back(v3);
                st_result[2].push_back(v4);
              }
              continue;
            }

            uint i6 = g_out_list[v5];
            for (const uintE *e6 = g_edge[i6]; i6 < g_out_list[v5 + 1]; ++i6) {
              const uint v6 = *e6;
              ++e6;
              const uint w6 = *e6;
              ++e6;
              if (!(info.dist[v6] & 4) || v6 == v1 || v6 == v2 || v6 == v3 ||
                  v6 == v4 || !cmpWeight(w5, w6))
                continue;
              else if (v6 == st) {
                if (cmpWeight(w6, w1)) {
                  st_result[3].push_back(v1);
                  st_result[3].push_back(v2);
                  st_result[3].push_back(v3);
                  st_result[3].push_back(v4);
                  st_result[3].push_back(v5);
                }
                continue;
              }
              const uint w7 = info.end_weight[v6];
              if (cmpWeight(w6, w7) && cmpWeight(w7, w1)) {
                st_result[4].push_back(v1);
                st_result[4].push_back(v2);
                st_result[4].push_back(v3);
                st_result[4].push_back(v4);
                st_result[4].push_back(v5);
                st_result[4].push_back(v6);
              }
            }
          }
        }
      }
    }
  }
}

void findCycleAt(const uint st, WorkerFindInfo &info) {
  backward3(st, info);
  if (info.visit_num < 2) return;
  forward7(st, info);
}

void subBFlower(const uint st, const uint left, const uint right,
                WorkerFindInfo &info) {
  for (int i = 0; i < info.visit_num; ++i) {
    info.dist[info.visit_list[i]] = 0;
  }
  info.visit_list[0] = st;
  info.visit_num = 1;
  info.dist[st] = 7;
  uint i1 = left;
  for (const uintE *e1 = g_reverse_edge[i1]; i1 < right; ++i1) {
    const uint v1 = *e1;
    ++e1;
    if (v1 <= st) break;
    const uint w1 = *e1;
    ++e1;
    info.dist[v1] = 7;  // 111
    info.end_weight[v1] = w1;
    info.visit_list[info.visit_num++] = v1;
    const uintE *e2 = g_reverse_edge[g_in_list[v1]];
    for (uint i2 = g_in_list[v1]; i2 < g_in_list[v1 + 1]; ++i2) {
      const uint v2 = *e2;
      ++e2;
      if (v2 <= st) break;
      const uint w2 = *e2;
      ++e2;
      if (!cmpWeight(w2, w1)) continue;
      info.dist[v2] |= 3;  // 011
      info.visit_list[info.visit_num++] = v2;
      const uintE *e3 = g_reverse_edge[g_in_list[v2]];
      for (uint i3 = g_in_list[v2]; i3 < g_in_list[v2 + 1]; ++i3) {
        const uint v3 = *e3;
        ++e3;
        if (v3 <= st) break;
        const uint w3 = *e3;  // HINT:位置不能错
        ++e3;
        if (v3 == v1) continue;
        if (!cmpWeight(w3, w2)) continue;
        info.dist[v3] |= 1;  // 001
        info.visit_list[info.visit_num++] = v3;
      }
    }
  }
}

void backwardFlower(const uint st) {
  uint left = g_in_list[st];
  uint right = g_in_list[st + 1];
  while (left < right &&
         g_reverse_edge[right - 1][0] < st) {  //不能换顺序 防止left==right==0
    --right;
  }

  uint step = (right - left) / kThreadNum;
  thread pool[kThreadNum - 1];
  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(subBFlower, st, left + step * i, left + step * (i + 1),
                     ref(w_find_info[i]));
  }
  subBFlower(st, left + step * (kThreadNum - 1), right,
             w_find_info[kThreadNum - 1]);
  for (auto &th : pool) th.join();
  // merge data
  // printf("st%d subed\n",st);
  auto &base_info = w_find_info[0];
  for (uint i = 1; i < kThreadNum; ++i) {
    for (uint j = 0; j < w_find_info[i].visit_num; ++j) {
      uint v = w_find_info[i].visit_list[j];
      if (base_info.dist[v] == 0)
        base_info.visit_list[base_info.visit_num++] = v;
      base_info.dist[v] |= w_find_info[i].dist[v];
      base_info.end_weight[v] |= w_find_info[i].end_weight[v];
    }
  }
}

void subFFlower(const uint st, const uint left, const uint right,
                WorkerFindInfo &info, vector<uint> (&st_result)[5]) {
  // st_result[4].reserve(400);
  uint i1 = left;
  for (const uintE *e1 = g_edge[i1]; i1 < right; ++i1) {
    const uint v1 = *e1;
    ++e1;
    const uint w1 = *e1;
    ++e1;
    if (v1 < st) continue;

    uint i2 = g_out_list[v1];
    for (const uintE *e2 = g_edge[i2]; i2 < g_out_list[v1 + 1]; ++i2) {
      const uint v2 = *e2;
      ++e2;
      const uint w2 = *e2;
      ++e2;
      if (v2 <= st || !cmpWeight(w1, w2)) continue;

      uint i3 = g_out_list[v2];
      for (const uintE *e3 = g_edge[i3]; i3 < g_out_list[v2 + 1]; ++i3) {
        const uint v3 = *e3;
        ++e3;
        const uint w3 = *e3;
        ++e3;
        if (v3 < st || v3 == v1 || !cmpWeight(w2, w3)) continue;
        if (v3 == st) {
          if (cmpWeight(w3, w1)) {
            // st_result[0].insert(st_result[0].end(), {st, v1, v2});
            // testnum++;
            // st_result[0].push_back(st);
            st_result[0].push_back(v1);
            st_result[0].push_back(v2);
          }
          continue;
        };
        uint i4 = g_out_list[v3];
        for (const uintE *e4 = g_edge[i4]; i4 < g_out_list[v3 + 1]; ++i4) {
          const uint v4 = *e4;
          ++e4;
          const uint w4 = *e4;
          ++e4;
          if (!(info.dist[v4] & 1) || v4 == v2 || v4 == v1 ||
              !cmpWeight(w3, w4))
            continue;
          else if (v4 == st) {
            if (cmpWeight(w4, w1)) {
              st_result[1].push_back(v1);
              st_result[1].push_back(v2);
              st_result[1].push_back(v3);
            }
            continue;
          }

          uint i5 = g_out_list[v4];
          for (const uintE *e5 = g_edge[i5]; i5 < g_out_list[v4 + 1]; ++i5) {
            const uint v5 = *e5;
            ++e5;
            const uint w5 = *e5;
            ++e5;
            if (!(info.dist[v5] & 2) || v5 == v1 || v5 == v2 || v5 == v3 ||
                !cmpWeight(w4, w5))
              continue;
            else if (v5 == st) {
              if (cmpWeight(w5, w1)) {
                st_result[2].push_back(v1);
                st_result[2].push_back(v2);
                st_result[2].push_back(v3);
                st_result[2].push_back(v4);
              }
              continue;
            }

            uint i6 = g_out_list[v5];
            for (const uintE *e6 = g_edge[i6]; i6 < g_out_list[v5 + 1]; ++i6) {
              const uint v6 = *e6;
              ++e6;
              const uint w6 = *e6;
              ++e6;
              if (!(info.dist[v6] & 4) || v6 == v1 || v6 == v2 || v6 == v3 ||
                  v6 == v4 || !cmpWeight(w5, w6))
                continue;
              else if (v6 == st) {
                if (cmpWeight(w6, w1)) {
                  st_result[3].push_back(v1);
                  st_result[3].push_back(v2);
                  st_result[3].push_back(v3);
                  st_result[3].push_back(v4);
                  st_result[3].push_back(v5);
                }
                continue;
              }
              const uint w7 = info.end_weight[v6];
              if (cmpWeight(w6, w7) && cmpWeight(w7, w1)) {
                st_result[4].push_back(v1);
                st_result[4].push_back(v2);
                st_result[4].push_back(v3);
                st_result[4].push_back(v4);
                st_result[4].push_back(v5);
                st_result[4].push_back(v6);
              }
            }
          }
        }
      }
    }
  }
}

void forwardFlower(const uint st, WorkerFindInfo &data) {
  uint left = g_out_list[st];
  uint right = g_out_list[st + 1];
  while (g_edge[left][0] < st && left < right) {
    ++left;
  }
  uint step = (right - left) / kThreadNum;
  thread pool[kThreadNum - 1];
  vector<uint> st_result[kThreadNum][5];
  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(subFFlower, st, left + step * i, left + step * (i + 1),
                     ref(data), ref(st_result[i]));
  }

  subFFlower(st, left + step * (kThreadNum - 1), right, data,
             st_result[kThreadNum - 1]);

  for (auto &th : pool) th.join();

  for (uint i = 0; i < kThreadNum; ++i) {
    for (uint j = 0; j < 5; ++j) {
      g_result[st].path[j].insert(g_result[st].path[j].end(),
                                  st_result[i][j].begin(),
                                  st_result[i][j].end());
    }
  }
}

void findCycleFlower(const uint st) {
  backwardFlower(st);
  // backward3(st,w_find_info[0]);
  forwardFlower(st, w_find_info[0]);
}

void findCycleBalance(uint pid) {
  uint next_node;
  // uint avg_deg_thresh=(g_edge_num/g_node_num+1)*g_node_num/50000;
  uint avg_deg_thresh = (g_edge_num / g_node_num + 1) * g_node_num / 50000;
  uint deg;
  while (true) {
    while (g_find_lock.test_and_set()) {
    }
    next_node = g_find_num < g_node_num ? g_find_num++ : g_node_num;
    g_find_lock.clear();
    if (next_node >= g_node_num) break;

    deg = (g_out_list[next_node + 1] - g_out_list[next_node]) *
          (g_node_num - next_node) / g_node_num;
    if (deg > avg_deg_thresh) {
      while (g_find_lock.test_and_set()) {
      };
      g_flower_point.push_back(next_node);
      g_find_lock.clear();
    } else
      findCycleAt(next_node, w_find_info[pid]);
  }
}

void findCycle() {
  g_result = vector<Result>(g_node_num);

  thread pool[kThreadNum - 1];
  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(findCycleBalance, i);
  }
  findCycleBalance(kThreadNum - 1);
  for (auto &th : pool) th.join();

  cout << g_flower_point.size() << endl;
  for (auto &st : g_flower_point) {
    findCycleFlower(st);
  }
}

/**写答案*/
/**异步IO，3个线程prepare，一个fwrite，memcpy根据12/8/4字节对齐可大幅度加速*/
void assignWriteNode(uint all_answer_len) {
  uint seted_buffer = 0;
  uint cnt = 0;
  uint visit_times = 0;

  uint step = all_answer_len / kBufferNum + 1;
  // g_buffer_step = step * kMaxIdStrLen;
  uint next_bar = step;
  for (uint i = 0; i < 5; ++i) {
    for (uint j = 0; j < g_node_num; ++j) {
      cnt += g_result[j]
                 .path[i]
                 .size();  // TODO 这里*t+3/t+2会更好 但是prepare复杂度变高
      visit_times++;
      if (cnt > next_bar) {
        w_write_info[seted_buffer + 1].start_idx = visit_times;
        w_write_info[seted_buffer].end_idx = visit_times;
        next_bar += step;
        seted_buffer++;
      }
    }
  }
  w_write_info[seted_buffer].end_idx = visit_times;
}

char *writePrePare() {
  uint R = 0;
  uint R2 = 0;
  for (uint t = 0; t < 5; t++) {
    uint temp_R = 0;
    for (uint i = 0; i < g_node_num; ++i) {
      temp_R += g_result[i].path[t].size();
    }
    R += temp_R / (t + 2);
    R2 += temp_R;
  }
  result_num = R;
  assignWriteNode(R2);
  // int *x;
  char *temp_c = new char[kMaxIdStrLen];
  sprintf(temp_c, "%d\n", R);

#ifdef LOCAL
  printf("%d cycles, \n", R);
#endif
  return temp_c;
}

void prepareBuffer(int pid) {
#ifdef LOCAL
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t3 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif

  uint job;
  while (true) {
    while (g_buffer_lock.test_and_set()) {
    }
    job = g_prepared_buffer_num < kBufferNum ? g_prepared_buffer_num++
                                             : kBufferNum;
    g_buffer_lock.clear();
    if (job >= kBufferNum) break;
    while (job > g_next_buffer + 4) {
      usleep(100);
    }

    w_write_info[job].answer_buffer = new char
        [result_num * 10 * kMaxIdStrLen /
         kBufferNum];  //(char*)malloc(result_num*10*kMaxIdStrLen/kThreadNum);
    char *p = w_write_info[job].answer_buffer;
    int mod_res = 0;
    uint st = w_write_info[job].start_idx;
    uint ed = w_write_info[job].end_idx;
    uint st_i = st / g_node_num;
    uint st_j = st % g_node_num;
    uint visit_times = st;
    for (uint i = st_i; i < 5; ++i) {
      if (visit_times >= ed) break;
      uint mod_val = i + 2;
      for (uint j = st_j; j < g_node_num; ++j) {
        if (visit_times >= ed) break;
        const auto &prefix_str = g_node_str[j];
        for (auto &x : g_result[j].path[i]) {
          if (mod_res == 0) {
            memcpy(p, prefix_str.val, kMaxIdStrLen);
            p += prefix_str.len;
          }
          const auto &str = g_node_str[x];
          memcpy(p, str.val, kMaxIdStrLen);
          p += str.len;
          ++mod_res;
          if (mod_res == mod_val) {
            mod_res = 0;
            *(p - 1) = '\n';
          }
        }
        visit_times++;
      }
      st_j = 0;
    }
    w_write_info[job].answer_len = p - w_write_info[job].answer_buffer;
    w_write_info[job].is_finish = 1;
  }
#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("process %d prepare buffer time %fs   \n", pid, t4 - t3);
#endif
}

//异步IO
void writeAsyn(char *result_num_str) {
#ifdef LOCAL
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t3 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif

  FILE *fd = fopen(RESULT, "w");
  // fallocate(fd, 0, 0, all_answer_len);
  fwrite(result_num_str, strlen(result_num_str), 1, fd);
  g_next_buffer = 0;
  while (g_next_buffer < kBufferNum) {
    while (!w_write_info[g_next_buffer].is_finish) {
      usleep(10);
    }
    fwrite(w_write_info[g_next_buffer].answer_buffer,
           w_write_info[g_next_buffer].answer_len, 1, fd);
    g_next_buffer++;
  }
  fclose(fd);

#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("write data time %fs   \n", t4 - t3);
#endif
}

void writeResult() {
  char *result_sum_str = writePrePare();
  thread pool[kThreadNum - 1];
  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(prepareBuffer, i + 1);
  }
  //    if (kThreadNum==1)
  //        prepareBuffer(0);
  writeAsyn(result_sum_str);
  for (auto &th : pool) th.join();
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
  findCycle();
#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t3 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("find cycle time %fs\n", t3 - t2);
#endif

  writeResult();
}

int main() { solve(); }