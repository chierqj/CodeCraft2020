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
#define kMaxE 3000024              //最大边数
#define kMaxN 6000024              //最大点数
#define kThreadNum 4               //进程/线程数量
#define kMaxIdStrLen 12            //整数表示的ID有多长
const uint kCannotReachFlg = 2147123456;

string g_test_file;     //输入文件
string g_predict_file;  //输出文件

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

typedef uint
    Path;  //输入变成每个点5个vector，所有路径节点依次放进去，取的时候每隔3/4/5/6/7即为一行
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
vector<Edge> g_in_list[kMaxN];         //邻接表

vector<vector<Path>> g_result[5];
uint g_answer_len[5][kMaxN];

WorkerData g_worker_data[kThreadNum];  //每个线程的独立数据

/********************************************third
 * party************************************************/
void *memcpy_tiny(void *dst, const void *src, uint len) {
  // register
  unsigned char *dd = (unsigned char *)dst + len;
  const unsigned char *ss = (const unsigned char *)src + len;
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

/**
 * 找环的思路：
 * 基本和初赛思路一样，先反向搜索inlist记录搜到的点的距离dist，同时保存这些点到changed_list,这里的距离改为用三位二进制表示了
 * 然后正向搜环
 */
void findCycleAt(uint st, WorkerData &data) {
  // TODO w*3 w*5可能溢出  用longlong或者预先存的时候 判断>INTMAX/3 INTMAX/5
  // 若大于则直接置为INTMAX

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
  data.dist[st] = 7;
  uint temp_size[5]{};
  // TODO 提前计算size st~v2
  const ushort sum_len0 = g_node_str[st].len;
  for (uint i1 = g_out_list[st].first; i1 < g_out_list[st].second; ++i1) {
    const BriefEdge &e1 = g_brief_edge[i1];
    const uint &v1 = e1.v;
    if (v1 < st) continue;
    const ulong &w1 = e1.w;
    const ushort sum_len1 = sum_len0 + g_node_str[v1].len;
    for (uint i2 = g_out_list[v1].first; i2 < g_out_list[v1].second; ++i2) {
      const BriefEdge &e2 = g_brief_edge[i2];
      const uint &v2 = e2.v;
      const ulong &w2 = e2.w;
      if (v2 <= st || !cmpWeight(w1, w2)) continue;
      const ushort sum_len2 = sum_len1 + g_node_str[v2].len;
      for (uint i3 = g_out_list[v2].first; i3 < g_out_list[v2].second; ++i3) {
        const BriefEdge &e3 = g_brief_edge[i3];
        const uint &v3 = e3.v;
        const ulong &w3 = e3.w;
        if (v3 < st || v3 == v1 || !cmpWeight(w2, w3)) continue;
        if (v3 == st) {
          if (cmpWeight(w3, w1)) {
            // g_result[0][st].emplace_back(st, v1, v2);
            g_result[0][st].insert(g_result[0][st].end(), {st, v1, v2});
            temp_size[0] += (sum_len2);
          }
          continue;
        };
        const ushort sum_len3 = sum_len2 + g_node_str[v3].len;
        for (uint i4 = g_out_list[v3].first; i4 < g_out_list[v3].second; ++i4) {
          const BriefEdge &e4 = g_brief_edge[i4];
          const uint &v4 = e4.v;
          if (!(data.dist[v4] & 1)) continue;
          const ulong &w4 = e4.w;
          if (!cmpWeight(w3, w4))
            continue;
          else if (v4 == st) {
            if (cmpWeight(w4, w1)) {
              // g_result[1][st].emplace_back(st, v1, v2, v3);
              g_result[1][st].insert(g_result[1][st].end(), {st, v1, v2, v3});
              temp_size[1] += (sum_len3);
            }
            continue;
          } else if (v4 == v2 || v4 == v1)
            continue;
          const ushort len4 = g_node_str[v4].len;
          for (uint i5 = g_out_list[v4].first; i5 < g_out_list[v4].second;
               ++i5) {
            const BriefEdge &e5 = g_brief_edge[i5];
            const uint &v5 = e5.v;
            if (!(data.dist[v5] & 2)) continue;
            const ulong &w5 = e5.w;
            if (!cmpWeight(w4, w5))
              continue;
            else if (v5 == st) {
              if (cmpWeight(w5, w1)) {  // TODO 这里好像不用判断了

                // g_result[2][st].emplace_back(st, v1, v2, v3, v4);
                g_result[2][st].insert(g_result[2][st].end(),
                                       {st, v1, v2, v3, v4});
                temp_size[2] += (sum_len3 + len4);
              }
              continue;
            } else if (v5 == v1 || v5 == v2 || v5 == v3)
              continue;
            const ushort len5 = g_node_str[v5].len;
            for (uint i6 = g_out_list[v5].first; i6 < g_out_list[v5].second;
                 ++i6) {
              const BriefEdge &e6 = g_brief_edge[i6];
              const uint &v6 = e6.v;
              const ulong &w6 = e6.w;
              if (!(data.dist[v6] & 4)) continue;
              if (!cmpWeight(w5, w6))
                continue;
              else if (v6 == st) {
                if (cmpWeight(w6, w1)) {
                  // g_result[3][st].emplace_back(st, v1, v2, v3, v4, v5);
                  g_result[3][st].insert(g_result[3][st].end(),
                                         {st, v1, v2, v3, v4, v5});
                  temp_size[3] += (sum_len3 + len4 + len5);
                }
                continue;
              } else if (v6 == v1 || v6 == v2 || v6 == v3 || v6 == v4)
                continue;
              const ulong &w7 = data.end_weight[v6];
              if (cmpWeight(w6, w7) && cmpWeight(w7, w1)) {
                // g_result[4][st].emplace_back(st, v1, v2, v3, v4, v5, v6);
                g_result[4][st].insert(g_result[4][st].end(),
                                       {st, v1, v2, v3, v4, v5, v6});
                temp_size[4] += (sum_len3 + len4 + len5 + g_node_str[v6].len);
              }
            }
          }
        }
      }
    }
  }
  for (int i = 0; i < data.changed_num; ++i) {
    data.dist[data.changed_list[i]] = 0;
  }
  data.dist[st] = 0;
  data.changed_num = 0;
  g_answer_len[0][st] = temp_size[0];
  g_answer_len[1][st] = temp_size[1];
  g_answer_len[2][st] = temp_size[2];
  g_answer_len[3][st] = temp_size[3];
  g_answer_len[4][st] = temp_size[4];
}

void findCycleSkip(uint pid) {
  for (uint i = pid; i < g_node_num; i += kThreadNum) {
    if (g_in_list[i].empty() || !g_out_list[i].second) continue;
    // printf("find%d\n",i);
    findCycleAt(i, g_worker_data[pid]);
  }
}

vector<uint> g_start_idx_list;
vector<uint> g_start_offset_list;
vector<uint> g_end_idx_list;

void setProcessWriteNodeNum(uint all_answer_len) {
  g_start_idx_list = {0};
  g_start_offset_list = {0};
  g_end_idx_list = {};

  uint step = all_answer_len / kThreadNum + 1;
  uint next_bar = step;
  uint cnt = 0;
  uint visit_times = 0;
  uint offset = 0;

  for (uint i = 0; i < 5; ++i) {
    for (uint j = 0; j < g_node_num; ++j) {
      offset += g_answer_len[i][j];
      cnt += g_result[i][j].size();
      visit_times++;
      if (cnt > next_bar) {
        g_start_idx_list.push_back(visit_times);
        g_end_idx_list.push_back(visit_times);
        g_start_offset_list.push_back(offset);
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

void writePrePare() {
  uint R = 0;
  uint R2 = 0;
  for (uint t = 0; t < 5; t++) {
    uint temp_R = 0;
    for (auto &i : g_result[t]) {
      temp_R += i.size();
    }
    R += temp_R / (t + 3);
    R2 += temp_R;
  }

  // int *x;
  char *temp_c = new char[20];

  ulong all_answer_len = 0;

  for (uint i = 0; i < 5; ++i) {
    for (uint j = 0; j < g_node_num; ++j) {
      all_answer_len += g_answer_len[i][j];
    }
  }
  setProcessWriteNodeNum(R2);

  sprintf(temp_c, "%d\n", R);
  all_answer_len += strlen(temp_c);
#ifdef LOCAL
  printf("%d cycles, %d characters\n", R, all_answer_len);
#endif

  int fd = open(g_predict_file.c_str(), O_RDWR | O_CREAT, 0666);

  fallocate(fd, 0, 0, all_answer_len);
  // ftruncate(fd,all_answer_len);

  char *answer_mmap = (char *)mmap(NULL, all_answer_len, PROT_READ | PROT_WRITE,
                                   MAP_SHARED, fd, 0);

  lseek(fd, all_answer_len - 1, SEEK_SET);
  write(fd, "\0", 1);
  close(fd);
  char *p = answer_mmap;
  for (int i = 0; i < strlen(temp_c); ++i) {
    *p = temp_c[i];
    p++;
  }
  g_path_mmap_begin = p;
}

// TODO： 分块写入
void writeResult(int pid) {
  char *p = g_path_mmap_begin + g_start_offset_list[pid];
  int mod_res = 0;
  uint st = g_start_idx_list[pid];
  uint ed = g_end_idx_list[pid];
  uint visit_times = 0;
  for (uint i = 0; i < 5; ++i) {
    for (uint j = 0; j < g_node_num; ++j) {
      if (visit_times >= st) {
        for (auto &x : g_result[i][j]) {
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
}

void findCycleTheading(uint pid) {
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
  for (auto &i : g_result) i = vector<vector<Path>>(g_node_num);

  thread pool[kThreadNum - 1];
  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(findCycleTheading, i);
  }
  findCycleTheading(kThreadNum - 1);
  for (auto &th : pool) th.join();

  exit(0);
}

int main() {
  g_test_file = "/data/test_data.txt";
  g_predict_file = "/projects/student/result.txt";
#ifdef LOCAL
  g_test_file = "testdata/standard.txt";
  g_test_file = "testdata/test_data_massive.txt";
  g_predict_file = "data/result.txt";
#endif
  solve();
}