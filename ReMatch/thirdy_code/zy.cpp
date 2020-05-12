#include <iostream>
#include <utility>
#include <vector>
//#include <sstream>
#include <fstream>
//#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <map>
#include <string>
//#include <stack>
#include <algorithm>
#include <cstdio>
#include <unordered_map>
//#include<unordered_set>
#include <fcntl.h>
#include <unistd.h>
//#include <sys/mman.h>
#include <thread>
//#include<fcntl.h>
//#include <mutex>
//#include <condition_variable>
#include <sys/time.h>
#include <sys/types.h>

#include <atomic>
#include <random>

using namespace std;
#define uint uint32_t
#define ulong uint32_t
#define ushort uint8_t

/********************************************third
 * party************************************************/

template <class T>
class Vec {
  uint elem_cnt;
  uint elem_capacity;
  T *elem;

 public:
  /** Default constructor */
  explicit Vec<T>(size_t s)
      : elem_cnt(0),
        elem_capacity(s),
        elem{static_cast<T *>(::operator new(sizeof(T) * elem_capacity))} {};

  explicit Vec<T>() : elem_cnt(0), elem_capacity(0), elem(nullptr){};

  Vec<T>(const Vec<T> &other)
      : elem_cnt(other.elem_cnt),
        elem_capacity(other.elem_capacity),
        elem(new T[elem_cnt]) {
    for (size_t i = 0; i < elem_cnt; ++i) elem[i] = other.elem[i];
  }

  Vec<T> &operator=(const Vec<T> &other) {
    for (size_t i = 0; i < elem_cnt; ++i) elem[i] = other.elem[i];
    elem_cnt = other.elem_cnt;
    elem_capacity = other.elem_capacity;
    return *this;
  }

  T *begin() { return elem; }

  T *begin() const { return elem; }

  T *end() { return (elem + elem_cnt); }

  T *end() const { return (elem + elem_cnt); }

  T &operator[](size_t n) { return elem[n]; }

  T &operator[](size_t n) const { return elem[n]; }

  size_t size() const { return elem_cnt; }

  /*    size_t capacity() const {
          return elem_capacity;
      }


      void push_back(T &value) {
          if (elem_cnt >= elem_capacity) {
              reserve(elem_capacity * 2+8);
          }
          elem[elem_cnt++] = value;
      }
      void push_back(const T &value) {
          if (elem_cnt >= elem_capacity) {
              reserve(elem_capacity * 2+8);
          }
          elem[elem_cnt++] = value;
      }

      void push_list(const T &v0,const T &v1,const T &v2){
          if (elem_cnt+3 >= elem_capacity) {
              reserve(elem_capacity * 2+8);
          }
          elem[elem_cnt++] = v0;
          elem[elem_cnt++] = v1;
          elem[elem_cnt++] = v2;
      }
      void push_list(const T &v0,const T &v1,const T &v2,const T &v3){
          if (elem_cnt+4 >= elem_capacity) {
              reserve(elem_capacity * 2+8);
          }
          elem[elem_cnt++] = v0;
          elem[elem_cnt++] = v1;
          elem[elem_cnt++] = v2;
          elem[elem_cnt++] = v3;
      }
      void push_list(const T &v0,const T &v1,const T &v2,const T &v3,const T
     &v4){ if (elem_cnt+5 >= elem_capacity) { reserve(elem_capacity * 2+8);
          }
          elem[elem_cnt++] = v0;
          elem[elem_cnt++] = v1;
          elem[elem_cnt++] = v2;
          elem[elem_cnt++] = v3;
          elem[elem_cnt++] = v4;
      }
      void push_list(const T &v0,const T &v1,const T &v2,const T &v3,const T
     &v4,const T &v5){ if (elem_cnt+6 >= elem_capacity) { reserve(elem_capacity
     * 2+8);
          }
          elem[elem_cnt++] = v0;
          elem[elem_cnt++] = v1;
          elem[elem_cnt++] = v2;
          elem[elem_cnt++] = v3;
          elem[elem_cnt++] = v4;
          elem[elem_cnt++] = v5;
      }
      void push_list(const T &v0,const T &v1,const T &v2,const T &v3,const T
     &v4,const T &v5,const T &v6){ if (elem_cnt+7 >= elem_capacity) {
              reserve(elem_capacity * 2+8);
          }
          T *p=elem+elem_cnt;
          *(p++)=v0;
          *(p++)=v1;
          *(p++)=v2;
          *(p++)=v3;
          *(p++)=v4;
          *(p++)=v5;
          *(p++)=v6;
          elem_cnt+=7;
      }*/

  template <typename... Args>
  void emplace_back(Args &&... args) noexcept {
    if (elem_cnt == elem_capacity) {
      reserve(elem_capacity * 2 + 8);
    }
    new (&elem[elem_cnt++]) T(std::forward<Args>(args)...);
  }

  void reserve(const size_t &n) {
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

  bool empty() { return elem_cnt == 0; }
};

/********************************************local********************************************************/

const uint kMaxResult = 22000128;     //最大结果数
const uint kMaxE = 3000024;           //最大边数
const uint kMaxN = 6000024;           //最大点数
const uint kThreadNum = 4;            //进程/线程数量
const uint kMaxIdStrLen = 11;         //整数表示的ID有多长
const uint kUselessFlg = 2147123456;  //标记没有出边的点

const uint kBufferNum = 128;

string g_test_file;     //输入文件
string g_predict_file;  //输出文件

typedef pair<uint, ulong>
    Edge;  // in_list邻接表的边的数据结构  u-v-w  in_list[v]=[{u1,w1},{u2,w2}]
/*typedef struct BriefEdge {
    uint first;
    ulong second;

    BriefEdge() = default;
} BriefEdge;*/
typedef pair<uint, ulong> BriefEdge;
typedef struct RawEdge {
  uint u;
  uint v;
  ulong w;

  RawEdge() = default;
  // RawEdge(uint x, uint y, ulong weight) : u(x), v(y), w(weight) {};
} RawEdge;  //原始输入的边的数据结构

typedef struct NodeIdstr {
  ushort len;
  char val[kMaxIdStrLen];
  NodeIdstr() = default;
  // NodeIdstr(const char *p, uint l) : val(p), len(l) {};
} Ids;  //存储节点id的数据结构

// typedef uint Path;
// //输入变成每个点5个vector，所有路径节点依次放进去，取的时候每隔3/4/5/6/7即为一行
// TODO 存的时候可以不存自己  3 4 5 6 7-> 2 3 4 5 6
typedef struct Path {
  vector<uint> path[5];
  // uint len[5];
} Cycle;
/*typedef struct CirclePath {
    uint path[7];//不初始化 初始化浪费时间

    CirclePath() = default;

    CirclePath(uint a, uint b, uint c) {
        path[0] = a;
        path[1] = b;
        path[2] = c;
    }

    CirclePath(uint a, uint b, uint c, uint d) {
        path[0] = a;
        path[1] = b;
        path[2] = c;
        path[3] = d;
    }

    CirclePath(uint a, uint b, uint c, uint d, uint e) {
        path[0] = a;
        path[1] = b;
        path[2] = c;
        path[3] = d;
        path[4] = e;
    }

    CirclePath(uint a, uint b, uint c, uint d, uint e, uint f) {
        path[0] = a;
        path[1] = b;
        path[2] = c;
        path[3] = d;
        path[4] = e;
        path[5] = f;
    }

    CirclePath(uint a, uint b, uint c, uint d, uint e, uint f, uint g) {
        path[0] = a;
        path[1] = b;
        path[2] = c;
        path[3] = d;
        path[4] = e;
        path[5] = f;
        path[6] = g;
    }
} Path;//存储结果环的数据结构*/
typedef struct WorkerData {
  uint changed_num = 0;  //记录反向搜索改了多少点的dist
  ulong end_weight
      [kMaxN];  //记录反向搜索第一层的点，用于比较最后一条边和第一条边的weight是否合法
  uint changed_list[kMaxN];  //改了哪些点的dist
  char dist[kMaxN];          //反向搜索得到的点离st的距离

} WorkerData;  //存储每个线程所使用的独立的数据结构

/*********************************load*/
uint g_edge_num = 0;      //边的数量
uint g_node_num = 0;      //点的数量
char *g_path_mmap_begin;  //输出文件的mmap指针

uint g_worker_node_num[kThreadNum]{};  //每个线程unique的点数
unordered_map<uint, uint>
    g_worker_id_map[kThreadNum];  //每个线程的idmap  id-index
uint g_worker_map_offset
    [kThreadNum]{};  //用于拼接各个线程的unique的点，计算好偏移量后各个线程可以并行处理
uint g_worker_edge_num[kThreadNum]{};  //每个线程的边数
RawEdge g_worker_raw_edge[kThreadNum][kMaxE];

uint g_part_worker_edge_num[kThreadNum][kThreadNum]{{}};
RawEdge g_part_worker_raw_edge[kThreadNum][kThreadNum]
                              [kMaxE];  // 4线程 每个线程存4组edge uid%4

uint g_unordered_node_id[kMaxN];  //原始输入的点的id值，只用于排序
Ids g_unordered_node_str[kMaxN];  // id的str和strlen
Ids g_node_str[kMaxN];            //排序过后的str和strlen
RawEdge g_raw_edge[kMaxE];        //原始输入的边
/***********************findcycle*/
// pair<uint, uint> g_out_list[kMaxN]{};//前向星存储
uint g_in_list[kMaxN]{};
uint g_out_list[kMaxN]{};

BriefEdge g_brief_edge[kMaxE];   //不存u
Vec<Edge> g_in_list_old[kMaxN];  //邻接表

BriefEdge g_brief_reverse_edge[kMaxE];

// pair<uint,uint>g_in_list[kMaxN]{};

vector<Cycle> g_result;
// uint g_answer_len[5][kMaxN];
WorkerData g_worker_data[kThreadNum];  //每个线程的独立数据

vector<uint> g_start_idx_list;
vector<uint> g_end_idx_list;  // TODO  换成数组
// uint g_buffer_step;
// char *g_result_num_char;
char g_answer_buffer[kBufferNum][kMaxResult / kBufferNum * kMaxIdStrLen * 7];
uint g_answer_len[kBufferNum];
bool g_prepare_buffer_finish[kBufferNum];
// std::atomic<uint> value(0);
atomic_flag g_write_lock = ATOMIC_FLAG_INIT;
uint g_prepared_buffer_num = 0;
atomic_flag g_find_lock = ATOMIC_FLAG_INIT;
uint g_processed_num = 0;
// TODO 预处理出边1的点
// TODO 强连通分量
/******************************************functions*********************************************************************/
/** 一个新的逻辑，如果一个点只出现在y中而不出现在x中，那它必然不构成任何环，因此我们只需要为rawedge中的u建立映射即可
 * 又由于边是排序过的，因此仅当遍历中发现不一样的u时才尝试插入到idmap
 *
 * TODO 两种str策略，一种是加在rawedge里，另一种是只存intid，后面再反解析成str
 * */
template <typename T>
inline void parseInteger(char *p, const char *q, T &x) {
  x = *p++ - '0';
  while (p != q) {
    x = (x << 3) + (x << 1) + (*p - '0');
    ++p;
  }
}

/*void addEdge(uint x, uint y, ulong w) {
    uint pid = x % kThreadNum;
    RawEdge &e = g_worker_raw_edge[pid][g_worker_edge_num[pid]];
    e.u = x;
    e.v = y;
    e.w = w;
    g_worker_edge_num[pid]++;
}*/
void addEdge(uint x, uint y, ulong w, uint pid) {
  uint mod_pid = x % kThreadNum;
  RawEdge &e = g_part_worker_raw_edge[pid][mod_pid]
                                     [g_part_worker_edge_num[pid][mod_pid]];
  e.u = x;
  e.v = y;
  e.w = w;
  g_part_worker_edge_num[pid][mod_pid]++;
}

void analyzeBuffer(char *buffer, uint max_buffer_size, uint pid) {
  const char split = ',';
  // const char split_line = '\n';//TODO 用<'0'判断
  char *p = buffer;
  char *q = buffer;
  uint x, y, x_len, y_len;
  ulong w;
  const char *x_str, *y_str;
  while (q - buffer < max_buffer_size) {
    while (*q != split) q++;

    parseInteger(p, q, x);
    p = ++q;
    while (*q != split) q++;

    parseInteger(p, q, y);
    p = ++q;

    while (*q != '\r' && *q != '\n') q++;
    parseInteger(p, q, w);
    if (*q == '\r') ++q;
    p = ++q;
    addEdge(x, y, w, pid);
  }
}

bool cmpEdge(const RawEdge &a, const RawEdge &b) {
  if (a.u != b.u)
    return a.u < b.u;
  else
    return a.v < b.v;
}

// process pid=0   0-x   4-x  8-x
void sortRawEdge(uint pid) {
  RawEdge *p = g_worker_raw_edge[pid];
  g_worker_edge_num[pid] = 0;
  for (uint i = 0; i < kThreadNum; ++i) {
    memcpy(p, g_part_worker_raw_edge[i][pid],
           g_part_worker_edge_num[i][pid] * sizeof(RawEdge));
    p += g_part_worker_edge_num[i][pid];
    g_worker_edge_num[pid] += g_part_worker_edge_num[i][pid];
  }
  sort(g_worker_raw_edge[pid], g_worker_raw_edge[pid] + g_worker_edge_num[pid],
       cmpEdge);
}

//只对u建立，不对v建立map
void establishIdMap(uint pid) {
  uint last = -1;

  for (uint i = 0; i < g_worker_edge_num[pid]; ++i) {
    RawEdge &e = g_worker_raw_edge[pid][i];
    if (e.u != last) {
      uint &worker_node_num = g_worker_node_num[pid];
      last = e.u;
      g_worker_id_map[pid].insert(make_pair(last, worker_node_num++));
    }
  }
}

void intToStr(uint x, Ids &str) {
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

//对各个线程的id_map，将他们的value(index)给加上offset,offset是从前往后到它的id_map[i].size()的和
//这样 所有id_map里的value都是unique的，并且是连续的
// id_map[0] {0,0} {4,1}
// id_map[1] {1,0} {5,1}
//->  id_map[1]  {1,2}  {5,3}
// key/id 存下来 未排序 g_unordered_node_id {0,4,1,5}
// intToStr  int(id) -> char[12]
// find(5)    id_map[5%4][5]=3  g_unordered_node_id[3]=5
void mergeUniqueId(uint pid) {
  const uint offset = g_worker_map_offset[pid];
  for (auto &p : g_worker_id_map[pid]) {
    p.second += offset;
    g_unordered_node_id[p.second] = p.first;
    intToStr(p.first, g_unordered_node_str
                          [p.second]);  // TODO:排序好id后再直接映射到g_node_str
  }
}

// unordered_id[ 25 18 23]   %4=1 2 3
// argsort [1 2 0]
// arg_reflect[2 0 1]
// g_node_str 排序过后的 idstr   ["18","23","25"]
uint *idArgSort() {
  uint *arg_list = new uint[g_node_num];
  uint *arg_reflect = new uint[g_node_num];
  for (uint i = 0; i < g_node_num; ++i) arg_list[i] = i;

  std::sort(arg_list, arg_list + g_node_num, [](int pos1, int pos2) {
    return (g_unordered_node_id[pos1] < g_unordered_node_id[pos2]);
  });
  for (uint i = 0; i < g_node_num; ++i) {
    arg_reflect[arg_list[i]] = i;
    g_node_str[i] = g_unordered_node_str[arg_list[i]];
  }
  return arg_reflect;
}

// pid=2  raw_edge [{18 21}  {18,23}  {18 25} ]
// id_map[2]=[{18,1}]  id_map[1]=[{25,0}]  id_map[3]=[{23,2}]
// id->unordered_id.index arg_reflect[2 0 1]
// unordered_id.index->ordered_id.index 在node_str里的index   edge
// [{0,-1},{0,1},{0,2}  ]
void mapEdgeUV(const uint *arg_reflect, uint pid) {
  for (uint i = 0; i < g_worker_edge_num[pid]; ++i) {
    RawEdge &e = g_worker_raw_edge[pid][i];

    uint x_id = arg_reflect[g_worker_id_map[pid][e.u]];
    e.u = x_id;

    uint y_pid = e.v % kThreadNum;
    uint y_id = kUselessFlg;
    auto it = g_worker_id_map[y_pid].find(e.v);
    if (it != g_worker_id_map[y_pid].end()) {
      y_id = arg_reflect[it->second];
    }
    e.v = y_id;
  }
}

// TODO 优先队列
uint getNext(uint (&head)[kThreadNum]) {
  uint min_val = 2147483647;
  uint min_idx = kThreadNum;
  for (uint i = 0; i < kThreadNum; ++i) {
    if (g_worker_raw_edge[i][head[i]].u < min_val &&
        head[i] < g_worker_edge_num[i]) {
      min_val = g_worker_raw_edge[i][head[i]].u;
      min_idx = i;
    }
  }
  return min_idx;
}

void mergeRawEdge() {
  /*
      g_edge_num = 0;
      for (uint i = 0; i < kThreadNum; ++i) {
          for (uint j = 0; j < g_worker_edge_num[i]; ++j) {
              RawEdge *e_worker = &g_worker_raw_edge[i][j];
              if (e_worker->v == kUselessFlg) {
                  continue;
              }
              RawEdge *e_master = &g_raw_edge[g_edge_num];
              memcpy(e_master, e_worker, sizeof(RawEdge));
              g_edge_num++;
          }
      }
  */

  g_edge_num = 0;
  uint head[kThreadNum] = {0};
  uint next = getNext(head);
  uint last_u;
  RawEdge *e_master = g_raw_edge;
  while (next < kThreadNum) {
    RawEdge *e_worker = &g_worker_raw_edge[next][head[next]];
    last_u = e_worker->u;
    while (last_u == e_worker->u) {
      if (e_worker->v != kUselessFlg) {
        memcpy(e_master, e_worker, sizeof(RawEdge));
        ++e_master;
      }
      ++e_worker;
    }
    head[next] = e_worker - g_worker_raw_edge[next];
    next = getNext(head);
  }
  g_edge_num = e_master - g_raw_edge;
}

// in_list_old邻接表  in_list[1]={0}  in_list[2]={0,1}
// out_list前向星  E=[{0,1} {0,2} {0,3} {1,2} {1,5}]   out_list[0]={0,2}
// out_list[1]={3,4}
// TODO 分成连续的块做  可能会更快
/**HINT!!!!!!!!!!!
 * 要避免出现out_list[u]=0的情况，还要考虑最后几个node全部为0的情况  不然有bug*/
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
      g_in_list_old[e->v].emplace_back(e->u, e->w);
      BriefEdge &b = g_brief_edge[i];
      b.first = e->v;
      b.second = e->w;
    }
    ++e;
  }
  g_out_list[g_node_num] = g_edge_num;  //必须加 否则最后一项不对
}
/**HINT!!!!!!!! */
void buildGraphFixTail() {
  g_out_list[g_node_num] = g_edge_num;
  uint tail = g_node_num - 1;
  uint first_unzero = g_raw_edge[0].u;

  for (uint i = tail; i > first_unzero; --i) {
    if (g_out_list[i] == 0) g_out_list[i] = g_out_list[i + 1];
  }
}

/*void buildGraph(uint pid) {
    uint u_last=2147483647;
    RawEdge *e=g_raw_edge;
    for (uint i = 0; i < g_edge_num; ++i) {
        if (e->u!=u_last) {
            if (e->u % kThreadNum == pid) {
                g_out_list[e->u] = i;
            }
        }
        if (e->v % kThreadNum == pid) {
            g_in_list_old[e->v].emplace_back(e->u, e->w);
            BriefEdge &b = g_brief_edge[i];
            b.first = e->v;
            b.second = e->w;
        }
        ++e;
    }
}*/

bool cmpAdjReverse(const Edge &a, const Edge &b) { return a.first > b.first; }

//降序排序in_list  out_list升序的
void sortAdjList(uint st, uint ed) {
  for (int i = st; i < ed; ++i) {
    if (g_in_list_old[i].size() > 1)
      sort(g_in_list_old[i].begin(), g_in_list_old[i].end(), cmpAdjReverse);
  }
}

void buildReverseEdge(uint st, uint ed, uint offset) {
  uint idx = offset;
  for (int i = st; i < ed; ++i) {
    g_in_list[i] = idx;
    for (auto &e : g_in_list_old[i]) {
      g_brief_reverse_edge[idx] = move(e);
      ++idx;
    }
    // g_in_list[i].second=idx;
  }
  g_in_list[ed] = idx;
}

void prepareBuildReverseEdge(vector<uint> &st_list, vector<uint> &ed_list,
                             vector<uint> &offset_list) {
  for (uint i = 0; i < kThreadNum; ++i) {
    st_list[i] = i * (g_node_num / kThreadNum);
    ed_list[i] = (i + 1) * (g_node_num / kThreadNum);
  }
  ed_list[kThreadNum - 1] = g_node_num;
  uint offset = 0;
  offset_list[0] = 0;
  for (uint i = 0; i < kThreadNum; ++i) {
    offset_list[i] = offset;
    if (i == kThreadNum - 1) break;
    for (uint j = st_list[i]; j < ed_list[i]; ++j) {
      offset += g_in_list_old[j].size();
    }
  }
}

// TODO:多线程解析
void loadPipeline(char *buffer, uint max_buffer_size) {
#ifdef LOCAL
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif
  thread pool[kThreadNum - 1];
  uint step = max_buffer_size / kThreadNum;
  char *p = buffer;
  char *q = buffer + step;
  for (uint i = 0; i < kThreadNum - 1; ++i) {
    if (q - buffer >= max_buffer_size) {
      pool[i] = thread(analyzeBuffer, p, 0, i);
      continue;
    }
    while ((*q) != '\n') {
      q++;
    }
    pool[i] = thread(analyzeBuffer, p, q - p + 1, i);
    p = ++q;
    q = p + max_buffer_size;
  }
  analyzeBuffer(p, buffer + max_buffer_size - p, kThreadNum - 1);
  for (auto &th : pool) th.join();

#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("load data time stage1 %fs\n", t2 - t1);
#endif

  //所有worker边排序

  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(sortRawEdge, i + 1);
  }
  sortRawEdge(0);
  for (auto &th : pool) th.join();
#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t3 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("load data time stage2 %fs\n", t3 - t2);
#endif

  //所有worker建立各自的ID映射
  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(establishIdMap, i + 1);
  }
  establishIdMap(0);
  for (auto &th : pool) th.join();
#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("load data time stage3 %fs\n", t4 - t3);
#endif

  //计算worker映射的偏移量并计算总点数
  g_worker_map_offset[0] = 0;
  g_node_num = g_worker_node_num[0];
  for (uint i = 0; i < kThreadNum - 1; ++i) {
    g_worker_map_offset[i + 1] =
        g_worker_map_offset[i] + g_worker_id_map[i].size();
    g_node_num += g_worker_node_num[i + 1];
  }

  //建立unordered_node_id
  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(mergeUniqueId, i + 1);
  }
  mergeUniqueId(0);
  for (auto &th : pool) th.join();
#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t5 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("load data time stage4 %fs\n", t5 - t4);
#endif

  //点排序
  uint *arg_reflect = idArgSort();
#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t6 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("load data time stage5 %fs\n", t6 - t5);
#endif

  //对worker边进行映射
  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(mapEdgeUV, arg_reflect, i + 1);
  }
  mapEdgeUV(arg_reflect, 0);
  for (auto &th : pool) th.join();
#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t7 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("load data time stage6 %fs\n", t7 - t6);
#endif

  //合并worker的边到g_raw_edge
  //得到了处理完的所有的边 (u_index,v_index,w)  u相同的边是连续的,v有序
  // u_index和u的序相同

  mergeRawEdge();

#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t8 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("load data time stage7 %fs\n", t8 - t7);
#endif

  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(buildGraph, i + 1);
  }
  buildGraph(0);
  for (auto &th : pool) th.join();
  buildGraphFixTail();
#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t9 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("load data time stage8 %fs\n", t9 - t8);
#endif

  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(sortAdjList, i * g_node_num / kThreadNum,
                     (i + 1) * g_node_num / kThreadNum);
  }
  sortAdjList((kThreadNum - 1) * g_node_num / kThreadNum, g_node_num);
  for (auto &th : pool) th.join();
#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t10 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("load data time stage9 %fs\n", t10 - t9);
#endif

  vector<uint> st_list(kThreadNum, 0);
  vector<uint> ed_list(kThreadNum, 0);
  vector<uint> offset_list(kThreadNum, 0);
  prepareBuildReverseEdge(st_list, ed_list, offset_list);

  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(buildReverseEdge, st_list[i], ed_list[i], offset_list[i]);
  }
  buildReverseEdge(st_list[kThreadNum - 1], ed_list[kThreadNum - 1],
                   offset_list[kThreadNum - 1]);
  for (auto &th : pool) th.join();
#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t11 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("load data time stage10 %fs\n", t11 - t10);
#endif

  /*    ofstream fout;
      fout.open("data/new.txt");
      for(uint i=0;i<g_node_num;++i){
          fout<<g_out_list[i]<<endl;
      }
      fout.close();*/
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
 * 5.多线程建图，入边按邻接表存，出边按前向星存(由于边已经被sort过了，这时前向星中可以只存边的开始和结束序号)
 * 6.入边多线程sort一下
 * **/
void loadTestData() {
  ifstream fin(g_test_file.c_str(), std::ios::binary);
  if (!fin) {
    printf("Cannot open test file\n");
    exit(0);
  }

  uint max_buffer_size = fin.seekg(0, std::ios::end).tellg();
  char *file_buffer = new char[max_buffer_size];  // vector<char>(MAXS);
  fin.seekg(0, std::ios::beg).read(file_buffer, max_buffer_size);
  fin.close();
  loadPipeline(file_buffer, max_buffer_size);
}

/**
 * 复赛新增的条件，判断相邻两条边的权重是否合法
 *
 */
const ulong MULT3_MAX = 715827882;
const ulong MULT5_MAX = 429496729;

inline bool cmpWeight(const ulong &x, const ulong &y) {
  return (y > MULT5_MAX || x <= (y + (y << 2))) &&
         (x > MULT3_MAX || (x + (x << 1)) >= y);
}

/**
 * 找环的思路：3+7
 * 基本和初赛思路一样，先反向搜索inlist记录搜到的点的距离dist，同时保存这些点到changed_list,这里的距离改为用三位二进制表示了
 * 然后正向搜环
 * 多线程节点按模4分配，图均匀的情况下性能较好，后续可以考虑改成automic动态分配
 */
void findCycleAt(uint st, WorkerData &data) {
  // TODO w*3 w*5可能溢出  用longlong或者预先存的时候 判断>INTMAX/3 INTMAX/5
  // 若大于则直接置为INTMAX

  for (int i = 0; i < data.changed_num; ++i) {
    data.dist[data.changed_list[i]] = 0;
  }
  data.changed_list[0] = st;
  data.changed_num = 1;
  data.dist[st] = 7;

  uint i1 = g_in_list[st];
  for (const BriefEdge *e1 = &g_brief_reverse_edge[i1]; i1 < g_in_list[st + 1];
       ++i1, ++e1) {
    const uint v1 = e1->first;
    if (v1 < st) break;
    const ulong w1 = e1->second;
    data.dist[v1] = 7;  // 111
    data.end_weight[v1] = w1;
    data.changed_list[data.changed_num++] = v1;
    uint i2 = g_in_list[v1];
    for (const BriefEdge *e2 = &g_brief_reverse_edge[i2];
         i2 < g_in_list[v1 + 1]; ++i2, ++e2) {
      const uint v2 = e2->first;
      if (v2 <= st) break;
      const ulong w2 = e2->second;
      if (!cmpWeight(w2, w1)) continue;
      data.dist[v2] |= 3;  // 011
      data.changed_list[data.changed_num++] = v2;
      uint i3 = g_in_list[v2];
      for (const BriefEdge *e3 = &g_brief_reverse_edge[i3];
           i3 < g_in_list[v2 + 1]; ++i3, ++e3) {
        // const BriefEdge &e3 = g_brief_edge[i3];
        const uint v3 = e3->first;
        if (v3 <= st) break;
        const ulong w3 = e3->second;
        if (!cmpWeight(w3, w2)) continue;
        data.dist[v3] |= 1;  // 001
        data.changed_list[data.changed_num++] = v3;
      }
    }
  }

  // uint temp_size[5]{};
  // TODO 提前计算size st~v2
  // const ushort sum_len0 = g_node_str[st].len;
  vector<uint>(&st_result)[5] = g_result[st].path;
  // uint(&st_resut_len)[5] = g_result[st].len;
  i1 = g_out_list[st];
  for (const BriefEdge *e1 = &g_brief_edge[i1]; i1 < g_out_list[st + 1];
       ++i1, ++e1) {
    // const BriefEdge &e1 = g_brief_edge[i1];
    const uint v1 = e1->first;
    if (v1 < st) continue;
    const ulong w1 = e1->second;
    // const ushort sum_len1 = sum_len0 + g_node_str[v1].len;

    uint i2 = g_out_list[v1];
    for (const BriefEdge *e2 = &g_brief_edge[i2]; i2 < g_out_list[v1 + 1];
         ++i2, ++e2) {
      // const BriefEdge &e2 = g_brief_edge[i2];
      const uint v2 = e2->first;
      const ulong w2 = e2->second;
      if (v2 <= st || !cmpWeight(w1, w2)) continue;
      // const ushort sum_len2 = sum_len1 + g_node_str[v2].len;

      uint i3 = g_out_list[v2];
      for (const BriefEdge *e3 = &g_brief_edge[i3]; i3 < g_out_list[v2 + 1];
           ++i3, ++e3) {
        // const BriefEdge &e3 = g_brief_edge[i3];
        const uint v3 = e3->first;
        const ulong w3 = e3->second;
        if (v3 < st || v3 == v1 || !cmpWeight(w2, w3)) continue;
        if (v3 == st) {
          if (cmpWeight(w3, w1)) {
            st_result[0].insert(st_result[0].end(), {st, v1, v2});
            // st_result[0].push_list(st, v1, v2);
            // temp_size[0] += (sum_len2);
          }
          continue;
        };
        // const ushort sum_len3 = sum_len2 + g_node_str[v3].len;

        uint i4 = g_out_list[v3];
        for (const BriefEdge *e4 = &g_brief_edge[i4]; i4 < g_out_list[v3 + 1];
             ++i4, ++e4) {
          // const BriefEdge &e4 = g_brief_edge[i4];
          const uint v4 = e4->first;
          if (!(data.dist[v4] & 1)) continue;
          const ulong w4 = e4->second;
          if (!cmpWeight(w3, w4))
            continue;
          else if (v4 == st) {
            if (cmpWeight(w4, w1)) {
              // g_result[1][st].emplace_back(st, v1, v2, v3);
              // g_result[1][st].insert(g_result[1][st].end(), {st, v1, v2,
              // v3});
              st_result[1].insert(st_result[1].end(), {st, v1, v2, v3});
              // st_result[1].push_list(st, v1, v2, v3);
              // temp_size[1] += (sum_len3);
            }
            continue;
          } else if (v4 == v2 || v4 == v1)
            continue;
          // const ushort len4 = g_node_str[v4].len;

          uint i5 = g_out_list[v4];
          for (const BriefEdge *e5 = &g_brief_edge[i5]; i5 < g_out_list[v4 + 1];
               ++i5, ++e5) {
            // const BriefEdge &e5 = g_brief_edge[i5];
            const uint v5 = e5->first;
            if (!(data.dist[v5] & 2)) continue;
            const ulong w5 = e5->second;
            if (!cmpWeight(w4, w5))
              continue;
            else if (v5 == st) {
              if (cmpWeight(w5, w1)) {  // TODO  前向搜的时候 有weight判断
                                        // 有没有办法通过前面的让这里不用判断

                // g_result[2][st].emplace_back(st, v1, v2, v3, v4);
                // g_result[2][st].insert(g_result[2][st].end(), {st, v1, v2,
                // v3, v4});
                st_result[2].insert(st_result[2].end(), {st, v1, v2, v3, v4});
                // st_result[2].push_list(st, v1, v2, v3,v4);
                // temp_size[2] += (sum_len3 + len4);
              }
              continue;
            } else if (v5 == v1 || v5 == v2 || v5 == v3)
              continue;
            // const ushort len5 = g_node_str[v5].len;

            uint i6 = g_out_list[v5];
            for (const BriefEdge *e6 = &g_brief_edge[i6];
                 i6 < g_out_list[v5 + 1]; ++i6, ++e6) {
              // const BriefEdge &e6 = g_brief_edge[i6];
              const uint v6 = e6->first;
              const ulong w6 = e6->second;
              if (!(data.dist[v6] & 4)) continue;
              if (!cmpWeight(w5, w6))
                continue;
              else if (v6 == st) {
                if (cmpWeight(w6, w1)) {
                  // g_result[3][st].emplace_back(st, v1, v2, v3, v4, v5);
                  // g_result[3][st].insert(g_result[3][st].end(), {st, v1, v2,
                  // v3, v4, v5});
                  st_result[3].insert(st_result[3].end(),
                                      {st, v1, v2, v3, v4, v5});
                  // st_result[3].push_list(st, v1, v2, v3,v4,v5);
                  // temp_size[3] += (sum_len3 + len4 + len5);
                }
                continue;
              } else if (v6 == v1 || v6 == v2 || v6 == v3 || v6 == v4)
                continue;
              const ulong w7 = data.end_weight[v6];
              if (cmpWeight(w6, w7) && cmpWeight(w7, w1)) {
                st_result[4].insert(st_result[4].end(),
                                    {st, v1, v2, v3, v4, v5, v6});
                // st_result[4].push_list(st, v1, v2, v3,v4,v5,v6);
                // temp_size[4] += (sum_len3 + len4 + len5 +
                // g_node_str[v6].len);
              }
            }
          }
        }
      }
    }
  }
}

void findCycleSkip(uint pid) {
  for (uint i = pid; i < g_node_num; i += kThreadNum) {
    if (g_in_list_old[i].empty()) continue;
    // printf("find%d\n",i);
    findCycleAt(i, g_worker_data[pid]);
  }
}

void findCycleBalance(uint pid) {
  uint next_node;
  while (true) {
    while (g_find_lock.test_and_set()) {
    }
    next_node = g_processed_num < g_node_num ? g_processed_num++ : g_node_num;
    g_find_lock.clear();
    if (next_node >= g_node_num) break;
    findCycleAt(next_node, g_worker_data[pid]);
  }
}

void findCycleThreading(uint pid) {
#ifdef LOCAL
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif
  // findCycleSkip(pid);
  findCycleBalance(pid);
#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("process %d find cycle time %fs\n", pid, t2 - t1);
#endif
}

/**
 * 写入数据的思路：
 * 计算所有答案的总节点数量，平均分给各个线程。各个线程根据result准备好要写入文件的char
 * buffer
 * 由于可能存在分界点在一个节点内部的情况，用一个visit变量记录内层for的次数，用它的值作为分界线
 *
 * */

void assignWriteNode(uint all_answer_len) {
  g_start_idx_list = vector<uint>(kBufferNum, 0);
  g_end_idx_list = vector<uint>(kBufferNum, 0);
  uint seted_buffer = 0;
  uint cnt = 0;
  uint visit_times = 0;

  uint step = all_answer_len / kBufferNum + 1;
  // g_buffer_step = step * kMaxIdStrLen;
  uint next_bar = step;
  for (uint i = 0; i < 5; ++i) {
    for (uint j = 0; j < g_node_num; ++j) {
      cnt += g_result[j].path[i].size();
      visit_times++;
      if (cnt > next_bar) {
        g_start_idx_list[seted_buffer + 1] = visit_times;
        g_end_idx_list[seted_buffer] = visit_times;

        // cout << visit_times << " " << next_bar << " " << cnt << endl;

        next_bar += step;
        seted_buffer++;
      }
    }
  }
  g_end_idx_list[seted_buffer] = visit_times;

#ifdef LOCAL
  // for (int i = 0; i < kBufferNum; i++)
  // printf("%d ", g_start_idx_list[i]);
  // printf("%d  total:%d\n", g_end_idx_list[kThreadNum - 1], visit_times);
#endif
}

char *writePrePare() {
  uint R = 0;
  uint R2 = 0;
  for (uint t = 0; t < 5; t++) {
    uint temp_R = 0;
    for (uint i = 0; i < g_node_num; ++i) {
      temp_R += g_result[i].path[t].size();
    }
    R += temp_R / (t + 3);
    R2 += temp_R;
  }

  // int *x;
  char *temp_c = new char[20];
  assignWriteNode(R2);
  sprintf(temp_c, "%d\n", R);

#ifdef LOCAL
  printf("%d cycles, \n", R);
#endif
  return temp_c;
}

// TODO： 分块写入
char *memcpy_tiny(char *dst, const char *src, uint len) {
  char *dd = (char *)dst + len;
  const char *ss = (const char *)src + len;
  switch (len) {
    case 12:
      *((int *)(dd - 12)) = *((int *)(ss - 12));
    case 8:
      *((int *)(dd - 8)) = *((int *)(ss - 8));
    case 4:
      *((int *)(dd - 4)) = *((int *)(ss - 4));
      break;

    case 11:
      *((int *)(dd - 11)) = *((int *)(ss - 11));
    case 7:
      *((int *)(dd - 7)) = *((int *)(ss - 7));
      *((int *)(dd - 4)) = *((int *)(ss - 4));
      break;
    case 3:
      *((short *)(dd - 3)) = *((short *)(ss - 3));
      dd[-1] = ss[-1];
      break;

    case 10:
      *((int *)(dd - 10)) = *((int *)(ss - 10));
    case 6:
      *((int *)(dd - 6)) = *((int *)(ss - 6));
    case 2:
      *((short *)(dd - 2)) = *((short *)(ss - 2));
      break;

    case 9:
      *((int *)(dd - 9)) = *((int *)(ss - 9));
    case 5:
      *((int *)(dd - 5)) = *((int *)(ss - 5));
    case 1:
      dd[-1] = ss[-1];
      break;
    case 0:
    default:
      break;
  }
  return dd;
}

// TODO 按点存str-> 按边存  7->4
void prepareBuffer(int pid) {
#ifdef LOCAL
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t3 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif

  // g_worker_answer_buffer[pid] = new char[g_buffer_step];
  uint job;
  while (true) {
    while (g_write_lock.test_and_set()) {
    }
    job = g_prepared_buffer_num < kBufferNum ? g_prepared_buffer_num++
                                             : kBufferNum;
    g_write_lock.clear();
    if (job >= kBufferNum) break;
    char *p = g_answer_buffer[job];
    int mod_res = 0;
    uint st = g_start_idx_list[job];
    uint ed = g_end_idx_list[job];
    uint st_i = st / g_node_num;
    uint st_j = st % g_node_num;
    uint visit_times = st;
    for (uint i = st_i; i < 5; ++i) {
      if (visit_times >= ed) break;
      uint mod_val = i + 3;
      for (uint j = st_j; j < g_node_num; ++j) {
        if (visit_times >= ed) break;
        for (auto &x : g_result[j].path[i]) {
          const auto &str = g_node_str[x];
          p = memcpy_tiny(p, str.val, str.len);
          // p += str.len;
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
    g_answer_len[job] = p - g_answer_buffer[job];
    g_prepare_buffer_finish[job] = true;
  }
#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("process %d prepare buffer time %fs   \n", pid, t4 - t3);
#endif
}

//异步写入
void writeAsyn(char *result_num_str) {
#ifdef LOCAL
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t3 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif

  FILE *fd = fopen(g_predict_file.c_str(), "w");
  // fallocate(fd, 0, 0, all_answer_len);
  fwrite(result_num_str, strlen(result_num_str), 1, fd);
  uint next_buffer = 0;
  while (next_buffer < kBufferNum) {
    while (g_prepare_buffer_finish[next_buffer] == false) {
      usleep(10);
    }
    fwrite(g_answer_buffer[next_buffer], g_answer_len[next_buffer], 1, fd);
    next_buffer++;
  }
  fclose(fd);

#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("write data time %fs   \n", t4 - t3);
#endif
}

void solve() {
#ifdef LOCAL
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif
  loadTestData();
#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("load data time %fs\n", t4 - t1);
#endif

  // setProcessNodeNum();
  //    for (auto &i : g_result)
  //        i = vector<vector<Path>>(g_node_num);
  g_result = vector<Cycle>(g_node_num);
  // g_result.reserve(g_node_num);

  thread pool[kThreadNum - 1];
  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(findCycleThreading, i);
  }
  findCycleThreading(kThreadNum - 1);
  for (auto &th : pool) th.join();

#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif
  char *result_sum_str = writePrePare();
#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("write prepare time %fs\n", t2 - t1);
#endif

#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif
  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(prepareBuffer, i + 1);
  }
  // prepareBuffer(0);

  writeAsyn(result_sum_str);
  for (auto &th : pool) th.join();
#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("prepare buffer time %fs\n", t2 - t1);
#endif

  // writeBuffer(result_sum_str);
  exit(0);
}

int main() {
  g_test_file = "/data/test_data.txt";
  g_predict_file = "/projects/student/result.txt";
#ifdef LOCAL
  g_test_file = "../data/gen/test_data.txt";
  // g_test_file = "testdata/test_data_massive.txt";
  // g_test_file = "testdata/1004812/test_data.txt";
  // g_predict_file = "data/result.txt";
  g_predict_file = "/dev/shm/result.txt";
#endif
  solve();
}