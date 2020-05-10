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
#include <sys/mman.h>
#include <unistd.h>

#include <thread>
//#include<fcntl.h>
//#include <mutex>
//#include <condition_variable>
#include <sys/time.h>
#include <sys/types.h>

#include <random>

using namespace std;
#define uint uint32_t
#define ulong uint32_t
#define ushort uint8_t

const uint kMaxResult = 20000008;  //最大结果数
const uint kMaxE = 3000024;        //最大边数
const uint kMaxN = 6000024;        //最大点数
const uint kThreadNum = 4;         //进程/线程数量
const uint kMaxIdStrLen = 12;      //整数表示的ID有多长
const uint kCannotReachFlg = 2147123456;

string g_test_file;     //输入文件
string g_predict_file;  //输出文件

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

typedef pair<uint, ulong> Edge;  //邻接表的边的数据结构
typedef struct BriefEdge {
  uint v;
  ulong w;

  BriefEdge() = default;
} BriefEdge;
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
typedef struct Path {
  vector<uint> path[5];
  // uint len[5];
} Path;
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
  bool is_finish = false;

} WorkerData;  //存储每个线程所使用的独立的数据结构

uint g_edge_num = 0;      //边的数量
uint g_node_num = 0;      //点的数量
char *g_path_mmap_begin;  //输出文件的mmap指针

uint g_worker_node_num[kThreadNum]{};
unordered_map<uint, uint> g_worker_id_map[kThreadNum];
uint g_worker_map_offset[kThreadNum]{};
uint g_worker_edge_num[kThreadNum]{};
RawEdge g_worker_raw_edge[kThreadNum][kMaxE];

uint g_unordered_node_id[kMaxN];  //原始输入的点的id值，只用于排序
Ids g_unordered_node_str[kMaxN];  // id的str和strlen
Ids g_node_str[kMaxN];            //排序过后的str和strlen
RawEdge g_raw_edge[kMaxE];        //原始输入的边
pair<uint, uint> g_out_list[kMaxN]{};  //前向星存储
BriefEdge g_brief_edge[kMaxE];         //不存u
Vec<Edge> g_in_list[kMaxN];            //邻接表

vector<Path> g_result;
// uint g_answer_len[5][kMaxN];
WorkerData g_worker_data[kThreadNum];  //每个线程的独立数据

/********************************************third
 * party************************************************/
/*void *memcpy_tiny(void *dst, const void *src, uint len) {
    //register
    unsigned char *dd = (unsigned char *) dst + len;
    const unsigned char *ss = (const unsigned char *) src + len;
    switch (len) {
        case 12:
            *((int *) (dd - 12)) = *((int *) (ss - 12));
        case 8:
            *((int *) (dd - 8)) = *((int *) (ss - 8));
        case 4:
            *((int *) (dd - 4)) = *((int *) (ss - 4));
            break;
        case 11:
            *((int *) (dd - 11)) = *((int *) (ss - 11));
        case 7:
            *((int *) (dd - 7)) = *((int *) (ss - 7));
            break;
        case 3:
            *((short *) (dd - 3)) = *((short *) (ss - 3));
            dd[-1] = ss[-1];
            break;
        case 10:
            *((int *) (dd - 10)) = *((int *) (ss - 10));
        case 6:
            *((int *) (dd - 6)) = *((int *) (ss - 6));
        case 2:
            *((short *) (dd - 2)) = *((short *) (ss - 2));
            break;
        case 9:
            *((int *) (dd - 9)) = *((int *) (ss - 9));
        case 5:
            *((int *) (dd - 5)) = *((int *) (ss - 5));
        case 1:
            dd[-1] = ss[-1];
            break;
        case 0:
        default:
            break;
    }
    return dd;
}*/

/********************************************local********************************************************/

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

void addEdge(uint x, uint y, ulong w) {
  uint pid = x % kThreadNum;
  RawEdge &e = g_worker_raw_edge[pid][g_worker_edge_num[pid]];
  e.u = x;
  e.v = y;
  e.w = w;
  g_worker_edge_num[pid]++;
}

bool cmpEdge(const RawEdge &a, const RawEdge &b) {
  if (a.u != b.u)
    return a.u < b.u;
  else
    return a.v < b.v;
}

void sortRawEdge(uint pid) {
  sort(g_worker_raw_edge[pid], g_worker_raw_edge[pid] + g_worker_edge_num[pid],
       cmpEdge);
}

void buildIdMap(uint pid) {
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

void buildUnorderdNodeArray(uint pid) {
  const uint offset = g_worker_map_offset[pid];
  for (auto &p : g_worker_id_map[pid]) {
    p.second += offset;
    g_unordered_node_id[p.second] = p.first;
    intToStr(p.first, g_unordered_node_str
                          [p.second]);  // TODO:排序好id后再直接映射到g_node_str
  }
}

uint *idArgSort() {
  uint *arg_list = new uint[g_node_num];
  uint *arg_reflact = new uint[g_node_num];
  for (uint i = 0; i < g_node_num; ++i) arg_list[i] = i;

  std::sort(arg_list, arg_list + g_node_num, [](int pos1, int pos2) {
    return (g_unordered_node_id[pos1] < g_unordered_node_id[pos2]);
  });
  for (uint i = 0; i < g_node_num; ++i) {
    arg_reflact[arg_list[i]] = i;
    g_node_str[i] = g_unordered_node_str[arg_list[i]];
  }
  return arg_reflact;
}

void mapEdgeUV(const uint *arg_reflect, uint pid) {
  for (uint i = 0; i < g_worker_edge_num[pid]; ++i) {
    RawEdge &e = g_worker_raw_edge[pid][i];

    uint x_id = arg_reflect[g_worker_id_map[pid][e.u]];
    e.u = x_id;

    uint y_pid = e.v % kThreadNum;
    uint y_id = kCannotReachFlg;
    auto it = g_worker_id_map[y_pid].find(e.v);
    if (it != g_worker_id_map[y_pid].end()) {
      y_id = arg_reflect[it->second];
    }
    e.v = y_id;
  }
}

void buildGraph(uint pid) {
  for (uint i = 0; i < g_edge_num; ++i) {
    RawEdge &e = g_raw_edge[i];
    if (e.u % kThreadNum == pid) {
      uint x_id = e.u;

      if (g_out_list[x_id].second == 0) {
        g_out_list[x_id].first = i;
      }
      g_out_list[x_id].second = i + 1;
    }
    if (e.v % kThreadNum == pid) {
      g_in_list[e.v].emplace_back(e.u, e.w);
      BriefEdge &b = g_brief_edge[i];
      b.v = e.v;
      b.w = e.w;
    }
  }
}

bool cmpAdjReverse(const Edge &a, const Edge &b) { return a.first > b.first; }

void sortAdjList(uint st, uint ed) {
  for (int i = st; i < ed; ++i) {
    if (g_in_list[i].size() > 1)
      sort(g_in_list[i].begin(), g_in_list[i].end(), cmpAdjReverse);
  }
}

void analyzeBuff(char *buffer, uint max_buffer_size) {
#ifdef LOCAL
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif
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
    addEdge(x, y, w);
  }
#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("load data time stage1 %fs\n", t2 - t1);
#endif

  //所有worker边排序
  thread pool[kThreadNum - 1];
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
    pool[i] = thread(buildIdMap, i + 1);
  }
  buildIdMap(0);
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
    pool[i] = thread(buildUnorderdNodeArray, i + 1);
  }
  buildUnorderdNodeArray(0);
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
  g_edge_num = 0;
  for (uint i = 0; i < kThreadNum; ++i) {
    for (uint j = 0; j < g_worker_edge_num[i]; ++j) {
      RawEdge *e_worker = &g_worker_raw_edge[i][j];
      if (e_worker->v == kCannotReachFlg) {
        continue;
      }
      RawEdge *e_master = &g_raw_edge[g_edge_num];
      memcpy(e_master, e_worker, sizeof(RawEdge));
      g_edge_num++;
    }
  }

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
}

/**
 * 读数据的思路：
 * 1.扫一遍所有数据，建立id_index_map，把所有边和点的id存下来(unordered_node_id)，先不加到临街表里
 * 2.把所有点根据id再排一次序，得到一个映射，其能保证映射后点的index和id的顺序完全相同
 * 3.多线程把所有边根据映射后的index加入邻接表，后续可以将所有id的比较替换为直接对index比较
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
  analyzeBuff(file_buffer, max_buffer_size);
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

inline bool judge(const BriefEdge *e1, const BriefEdge *e2) {
  return (e2->w > MULT5_MAX || e1->w <= (e2->w + (e2->w << 2))) &&
         (e1->w > MULT3_MAX || (e1->w + (e1->w << 1)) >= e2->w);
}

/**
 * 找环的思路：
 * 基本和初赛思路一样，先反向搜索inlist记录搜到的点的距离dist，同时保存这些点到changed_list,这里的距离改为用三位二进制表示了
 * 然后正向搜环
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

  for (auto &e1 : g_in_list[st]) {
    const uint &v1 = e1.first;
    if (v1 < st) break;
    const ulong &w1 = e1.second;
    data.dist[v1] = 7;  // 111
    data.end_weight[v1] = w1;
    data.changed_list[data.changed_num++] = v1;
    for (auto &e2 : g_in_list[v1]) {
      const uint &v2 = e2.first;
      if (v2 <= st) break;
      const ulong &w2 = e2.second;
      if (!cmpWeight(w2, w1)) continue;
      data.dist[v2] |= 3;  // 011
      data.changed_list[data.changed_num++] = v2;
      for (auto &e3 : g_in_list[v2]) {
        const uint &v3 = e3.first;
        if (v3 <= st) break;
        const ulong &w3 = e3.second;
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
  uint i1 = g_out_list[st].first;
  for (const BriefEdge *e1 = &g_brief_edge[i1]; i1 < g_out_list[st].second;
       ++i1, ++e1) {
    // const BriefEdge &e1 = g_brief_edge[i1];
    const uint &v1 = e1->v;
    if (v1 < st) continue;
    const ulong &w1 = e1->w;
    // const ushort sum_len1 = sum_len0 + g_node_str[v1].len;

    uint i2 = g_out_list[v1].first;
    for (const BriefEdge *e2 = &g_brief_edge[i2]; i2 < g_out_list[v1].second;
         ++i2, ++e2) {
      // const BriefEdge &e2 = g_brief_edge[i2];
      const uint &v2 = e2->v;
      const ulong &w2 = e2->w;
      if (v2 <= st || !cmpWeight(w1, w2)) continue;
      // const ushort sum_len2 = sum_len1 + g_node_str[v2].len;

      uint i3 = g_out_list[v2].first;
      for (const BriefEdge *e3 = &g_brief_edge[i3]; i3 < g_out_list[v2].second;
           ++i3, ++e3) {
        // const BriefEdge &e3 = g_brief_edge[i3];
        const uint &v3 = e3->v;
        const ulong &w3 = e3->w;
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

        uint i4 = g_out_list[v3].first;
        for (const BriefEdge *e4 = &g_brief_edge[i4];
             i4 < g_out_list[v3].second; ++i4, ++e4) {
          // const BriefEdge &e4 = g_brief_edge[i4];
          const uint &v4 = e4->v;
          if (!(data.dist[v4] & 1)) continue;
          const ulong &w4 = e4->w;
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

          uint i5 = g_out_list[v4].first;
          for (const BriefEdge *e5 = &g_brief_edge[i5];
               i5 < g_out_list[v4].second; ++i5, ++e5) {
            // const BriefEdge &e5 = g_brief_edge[i5];
            const uint &v5 = e5->v;
            if (!(data.dist[v5] & 2)) continue;
            const ulong &w5 = e5->w;
            if (!cmpWeight(w4, w5))
              continue;
            else if (v5 == st) {
              if (cmpWeight(w5, w1)) {  // TODO 这里好像不用判断了

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

            uint i6 = g_out_list[v5].first;
            for (const BriefEdge *e6 = &g_brief_edge[i6];
                 i6 < g_out_list[v5].second; ++i6, ++e6) {
              // const BriefEdge &e6 = g_brief_edge[i6];
              const uint v6 = e6->v;
              const ulong w6 = e6->w;
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
                // g_result[4][st].emplace_back(st, v1, v2, v3, v4, v5, v6);
                // g_result[4][st].insert(g_result[4][st].end(), {st, v1, v2,
                // v3, v4, v5, v6});
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

  /*   uint sz0 = 0, sz1 = 0, sz2 = 0, sz3 = 0, sz4 = 0;
     const uint &len0 = g_node_str[st].len;
     auto &ret = g_result[st];
     const auto &cdr1 = g_out_list[st];
     const BriefEdge *e1 = &g_brief_edge[cdr1.first];
     for (uint it1 = cdr1.first; it1 != cdr1.second; ++it1, ++e1) {
         if (e1->v < st) continue;
         const uint &len1 = g_node_str[e1->v].len + len0;
         const auto &cdr2 = g_out_list[e1->v];
         const BriefEdge *e2 = &g_brief_edge[cdr2.first];
         for (uint it2 = cdr2.first; it2 != cdr2.second; ++it2, ++e2) {
             if (e2->v <= st || !judge(e1, e2)) continue;
             const uint &len = g_node_str[e2->v].len + len1;
             const auto &cdr3 = g_out_list[e2->v];
             const BriefEdge *e3 = &g_brief_edge[cdr3.first];
             for (uint it3 = cdr3.first; it3 != cdr3.second; ++it3, ++e3) {
                 if (e3->v < st || e3->v == e1->v || !judge(e2, e3)) {
                     continue;
                 } else if (e3->v == st) {
                     if (!judge(e3, e1)) continue;
                     ret.path[0].insert(ret.path[0].end(), {st, e1->v, e2->v});
                     sz0 += len;
                     continue;
                 }
                 const uint &len3 = g_node_str[e3->v].len;
                 const auto &cdr4 = g_out_list[e3->v];
                 const BriefEdge *e4 = &g_brief_edge[cdr4.first];
                 for (uint it4 = cdr4.first; it4 != cdr4.second; ++it4, ++e4) {
                     if (!(data.dist[e4->v] & 1) || e1->v == e4->v || e2->v ==
     e4->v) { continue; } else if (!judge(e3, e4)) { continue; } else if (e4->v
     == st) { if (!judge(e4, e1)) continue;
                         ret.path[1].insert(ret.path[1].end(), {st, e1->v,
     e2->v, e3->v}); sz1 += len + len3; continue;
                     }
                     const uint &len4 = g_node_str[e4->v].len;
                     const auto &cdr5 = g_out_list[e4->v];
                     const BriefEdge *e5 = &g_brief_edge[cdr5.first];
                     for (uint it5 = cdr5.first; it5 != cdr5.second; ++it5,
     ++e5) { if (!(data.dist[e5->v] & 2) || e1->v == e5->v || e2->v == e5->v ||
                             e3->v == e5->v) {
                             continue;
                         } else if (!judge(e4, e5)) {
                             continue;
                         } else if (e5->v == st) {
                             if (!judge(e5, e1)) continue;
                             ret.path[2].insert(ret.path[2].end(),
                                               {st, e1->v, e2->v, e3->v,
     e4->v}); sz2 += len + len3 + len4; continue;
                         }
                         const uint &len5 = g_node_str[e5->v].len;
                         const auto &cdr6 = g_out_list[e5->v];
                         const BriefEdge *e6 = &g_brief_edge[cdr6.first];
                         for (uint it6 = cdr6.first; it6 != cdr6.second; ++it6,
     ++e6) { if (!(data.dist[e6->v] & 4) || e1->v == e6->v || e2->v == e6->v ||
     e3->v == e6->v || e4->v == e6->v) { continue; } else if (!judge(e5, e6)) {
                                 continue;
                             } else if (e6->v == st) {
                                 if (!judge(e6, e1)) continue;
                                 ret.path[3].insert(ret.path[3].end(),
                                                   {st, e1->v, e2->v, e3->v,
     e4->v, e5->v}); sz3 += len + len3 + len4 + len5; continue;
                             }
                             const uint &w7 = data.end_weight[e6->v];
                             if (!cmpWeight(e6->w, w7) || !cmpWeight(w7, e1->w))
     { continue;
                             }
                             const uint &len6 = g_node_str[e6->v].len;
                             ret.path[4].insert(ret.path[4].end(),
                                               {st, e1->v, e2->v, e3->v, e4->v,
     e5->v, e6->v}); sz4 += len + len3 + len4 + len5 + len6;
                         }
                     }
                 }
             }
         }
     }

     ret.len[0] = sz0;
     ret.len[1] = sz1;
     ret.len[2] = sz2;
     ret.len[3] = sz3;
     ret.len[4] = sz4;*/
}

void findCycleSkip(uint pid) {
  for (uint i = pid; i < g_node_num; i += kThreadNum) {
    if (g_in_list[i].empty()) continue;
    // printf("find%d\n",i);
    findCycleAt(i, g_worker_data[pid]);
  }
}

vector<uint> g_start_idx_list;
vector<uint> g_end_idx_list;
uint g_buffer_step;

void setProcessWriteNodeNum(uint all_answer_len) {
  g_start_idx_list = {0};
  g_end_idx_list = {};

  uint step = all_answer_len / kThreadNum + 1;
  // g_buffer_step = step * kMaxIdStrLen;
  uint next_bar = step;
  uint cnt = 0;
  uint visit_times = 0;

  for (uint i = 0; i < 5; ++i) {
    for (uint j = 0; j < g_node_num; ++j) {
      cnt += g_result[j].path[i].size();
      visit_times++;
      if (cnt > next_bar) {
        g_start_idx_list.push_back(visit_times);
        g_end_idx_list.push_back(visit_times);
        cout << i << " " << next_bar << " " << cnt << endl;
        next_bar += step;
      }
    }
  }

  /*    while (g_start_idx_list.size() < kThreadNum) {
          g_start_idx_list.push_back(g_node_num);
          g_end_idx_list.push_back(g_node_num);
      }*/
  g_end_idx_list.push_back(visit_times);
#ifdef LOCAL
  for (int i = 0; i < kThreadNum; i++) printf("%d ", g_start_idx_list[i]);
  printf("%d  total:%d\n", g_end_idx_list[kThreadNum - 1], visit_times);
#endif
}

char *g_result_num_char;
char *g_worker_answer_buffer[kThreadNum];
uint g_worker_answer_len[kThreadNum];
void writePrePare() {
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

  setProcessWriteNodeNum(R2);

  sprintf(temp_c, "%d\n", R);
  g_result_num_char = temp_c;

#ifdef LOCAL
  printf("%d cycles, \n", R);
#endif

  // memset(answer_mmap, '\0', all_answer_len);
}

// TODO： 分块写入
void writeResult(int pid) {
  g_worker_answer_buffer[pid] = new char[g_buffer_step];
  char *p = g_worker_answer_buffer[pid];
  int mod_res = 0;
  uint st = g_start_idx_list[pid];
  uint ed = g_end_idx_list[pid];
  uint visit_times = 0;
  for (uint i = 0; i < 5; ++i) {
    for (uint j = 0; j < g_node_num; ++j) {
      if (visit_times >= st) {
        for (auto &x : g_result[j].path[i]) {
          const auto &str = g_node_str[x];
          memcpy(p, str.val, str.len);
          p += str.len;
          ++mod_res;
          if (mod_res == i + 3) {
            mod_res = 0;
            *(p - 1) = '\n';
          }
        }
      }
      visit_times++;
      if (visit_times >= ed) break;
    }
  }
  g_worker_answer_len[pid] = p - g_worker_answer_buffer[pid];
}

/*void writeResult(int pid) {
    char *p = g_path_mmap_begin;
    int mod_res = 0;
    for (uint i = 0; i < 5; ++i) {
        for (uint j = 0; j < g_node_num; ++j) {
            if (j % kThreadNum == pid) {
                for (auto &x:g_result[j].path[i]) {
                    const auto &str = g_node_str[x];
                    memcpy(p, str.val, str.len);
                    p += str.len;
                    ++mod_res;
                    if (mod_res == i + 3) {
                        mod_res = 0;
                        *(p - 1) = '\n';
                    }
                }
            } else {
                p += g_result[j].len[i];
            }
        }
    }

}*/

void findCycleThreading(uint pid) {
#ifdef LOCAL
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif

  findCycleSkip(pid);
#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("process %d find cycle time %fs\n", pid, t2 - t1);
#endif
  if (pid == 0) {
    while (true) {
      bool flg = true;
      for (uint i = 1; i < kThreadNum; i++) {
        flg = flg && g_worker_data[i].is_finish;
      }
      if (flg) break;
      usleep(10);
    }
#ifdef LOCAL
    gettimeofday(&tim, nullptr);
    t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif
    writePrePare();
#ifdef LOCAL
    gettimeofday(&tim, nullptr);
    t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
    printf("process %d write prepare time %fs\n", pid, t2 - t1);
#endif
  }
  g_worker_data[pid].is_finish = true;

  while (true) {
    bool flg = true;
    for (auto &i : g_worker_data) {
      flg = flg && i.is_finish;
    }
    if (flg) break;
    usleep(10);
  }
#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t3 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif

  writeResult(pid);

#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("process %d write data time %fs   wait%fs\n", pid, t4 - t3, t3 - t2);
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
  g_result = vector<Path>(g_node_num);
  // g_result.reserve(g_node_num);

  thread pool[kThreadNum - 1];
  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(findCycleThreading, i);
  }
  findCycleThreading(kThreadNum - 1);
  for (auto &th : pool) th.join();

  FILE *fd = fopen(g_predict_file.c_str(), "w+");
  // fallocate(fd, 0, 0, all_answer_len);
  fwrite(g_result_num_char, strlen(g_result_num_char), 1, fd);
  for (uint i = 0; i < kThreadNum; ++i) {
    fwrite(g_worker_answer_buffer[i], g_worker_answer_len[i], 1, fd);
  }

  exit(0);
}

int main() {
  g_test_file = "/data/test_data.txt";
  g_predict_file = "/projects/student/result.txt";
#ifdef LOCAL
  g_test_file = "testdata/standard.txt";
  g_test_file = "testdata/test_data_massive.txt";
  // g_test_file = "testdata/1004812/test_data.txt";
  g_predict_file = "data/result.txt";
  g_predict_file = "/tmp/result.txt";
#endif
  solve();
}