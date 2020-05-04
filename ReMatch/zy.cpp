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

#include <random>

using namespace std;
#define uint uint32_t
#define ulong uint32_t

string g_test_file;                //输入文件
string g_predict_file;             //输出文件
const uint kMaxResult = 20000008;  //最大结果数
const uint kMaxE = 3000008;        //最大边数
const uint kMaxN = 6000008;        //最大点数
const uint kThreadNum = 4;         //进程/线程数量
const uint kMaxIdStrLen = 11;      //整数表示的ID有多长

typedef pair<uint, ulong> Edge;  //邻接表的边的数据结构
// typedef uint Edge;
// struct Edge {
//    uint first;
//    uint eid;
//
//    Edge() = default;
//
//    Edge(uint x,  uint e) : first(x), eid(e) {};
//};

typedef struct RawEdge {
  const char *val;
  uint u;
  uint v;
  ulong w;
  char len;
  char half_len;  // TODO short 貌似略快

  RawEdge() = default;

  // RawEdge(uint x, uint y, ulong weight) : u(x), v(y), w(weight) {};
} RawEdge;  //原始输入的边的数据结构

typedef struct CirclePath {
  uint path[4];  //不初始化 初始化浪费时间

  CirclePath() = default;

  CirclePath(int a, int b) {
    path[0] = a;
    path[1] = b;
  }

  CirclePath(int a, int b, int c) {
    path[0] = a;
    path[1] = b;
    path[2] = c;
  }

  CirclePath(int a, int b, int c, int d) {
    path[0] = a;
    path[1] = b;
    path[2] = c;
    path[3] = d;
  }
} Path;  //存储结果环的数据结构

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

uint g_unordered_node_id[kMaxN];  //原始输入的点的id值，只用于排序

RawEdge g_raw_edge[kMaxE];             //原始输入的边
vector<Edge> g_in_list[kMaxN];         //邻接表
pair<uint, uint> g_out_list[kMaxN]{};  //前向星存储

vector<vector<Path>> g_result[5];
// vector<vector<uint>> g_result[5];
uint g_answer_len[5][kMaxN];

WorkerData g_worker_data[kThreadNum];  //每个线程的独立数据

/********************************************third
 * party************************************************/
void *memcpy_tiny(void *dst, const void *src, size_t len) {
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

uint *idArgSort() {
  uint *arg_list = new uint[g_node_num];
  uint *arg_reflact = new uint[g_node_num];
  for (uint i = 0; i < g_node_num; ++i) arg_list[i] = i;

  std::sort(arg_list, arg_list + g_node_num, [](int pos1, int pos2) {
    return (g_unordered_node_id[pos1] < g_unordered_node_id[pos2]);
  });
  for (uint i = 0; i < g_node_num; ++i) {
    arg_reflact[arg_list[i]] = i;
    // g_node_ids[i] = g_unordered_node_ids[arg_list[i]];
  }
  return arg_reflact;
}

void buildGraph(const uint *arg_reflect,
                unordered_map<uint, uint> &id_index_map, uint pid) {
  for (uint i = 0; i < g_edge_num; ++i) {
    RawEdge &e = g_raw_edge[i];
    if (e.u % kThreadNum == pid) {
      uint x_id = arg_reflect[id_index_map[e.u]];
      e.u = x_id;
      if (g_out_list[x_id].second == 0) {
        g_out_list[x_id].first = i;
      }
      g_out_list[x_id].second++;
    }
    if (e.v % kThreadNum == pid) {
      ulong w = e.w;
      uint x_id = arg_reflect[id_index_map[e.u]];
      uint y_id = arg_reflect[id_index_map[e.v]];
      g_in_list[y_id].emplace_back(x_id, w);
    }
  }
}

bool cmpEdge(const RawEdge &a, const RawEdge &b) {
  if (a.u != b.u)
    return a.u < b.u;
  else
    return a.v < b.v;
}

bool cmpAdjReverse(const Edge &a, const Edge &b) { return a.first > b.first; }

void sortAdjList(uint st, uint ed) {
  for (int i = st; i < ed; ++i) {
    if (g_in_list[i].size() > 1)
      sort(g_in_list[i].begin(), g_in_list[i].end(), cmpAdjReverse);
  }
}

template <typename T>
inline void parseInteger(char *p, const char *q, T &x) {
  x = *p++ - '0';
  while (p != q) {
    x = (x << 3) + (x << 1) + (*p - '0');
    ++p;
  }
}

void addEdge(const char *str, uint len, uint half_len, uint x, uint y, ulong w,
             unordered_map<uint, uint> &id_index_map) {
  auto it_x = id_index_map.find(x);
  auto it_y = id_index_map.find(y);
  auto it_ed = id_index_map.end();

  if (it_x == it_ed) {
    g_unordered_node_id[g_node_num] = x;
    id_index_map.insert(make_pair(x, g_node_num));
    ++g_node_num;
    it_ed = id_index_map.end();
  }

  if (it_y == it_ed) {
    g_unordered_node_id[g_node_num] = y;
    id_index_map.insert(make_pair(y, g_node_num));
    ++g_node_num;
  }
  RawEdge &e = g_raw_edge[g_edge_num];
  e.u = x;
  e.v = y;
  e.w = w;
  e.val = str;
  e.len = len;
  e.half_len = half_len;
  ++g_edge_num;
}

void analyzeBuff(char *buffer, uint max_buffer_size) {
  const char split = ',';
  const char split_line = '\n';  // TODO 线上数据为/r/n
  char *p = buffer;
  char *q = buffer;

  uint len;

  uint x, y, x_len, y_len;
  ulong w;
  const char *x_str, *y_str;
  unordered_map<uint, uint> id_index_map;
  id_index_map.reserve(1000000);
  id_index_map.rehash(1000000);

  while (q - buffer < max_buffer_size) {
    while (*q != split) q++;
    // x_str = p;
    x_str = p;
    x_len = q - p + 1;

    parseInteger(p, q, x);

    p = ++q;
    while (*q != split) q++;

    // y_str = p;
    // y_len = q - p + 1;
    parseInteger(p, q, y);

    p = ++q;

    len = q - x_str;

    while (*q != '\r' && *q != '\n') q++;
    parseInteger(p, q, w);
    if (*q == '\r') ++q;
    p = ++q;
    // addEdge(x_str, y_str, x_len, y_len, x, y, w, id_index_map);
    addEdge(x_str, len, x_len, x, y, w, id_index_map);
  }

  cout << g_edge_num << " " << id_index_map.size() << endl;

  sort(g_raw_edge, g_raw_edge + g_edge_num, cmpEdge);

  uint *arg_reflect = idArgSort();

  thread pool[kThreadNum - 1];

  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(buildGraph, arg_reflect, ref(id_index_map), i + 1);
  }
  buildGraph(arg_reflect, ref(id_index_map), 0);
  for (auto &th : pool) th.join();

  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(sortAdjList, i * g_node_num / kThreadNum,
                     (i + 1) * g_node_num / kThreadNum);
  }
  sortAdjList((kThreadNum - 1) * g_node_num / kThreadNum, g_node_num);
  for (auto &th : pool) th.join();
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
  // TODO 提前计算size st~v2
  for (uint i1 = g_out_list[st].first;
       i1 < g_out_list[st].first + g_out_list[st].second; ++i1) {
    const RawEdge &e1 = g_raw_edge[i1];
    const uint &v1 = e1.v;
    const ulong &w1 = e1.w;
    if (v1 < st) continue;
    for (uint i2 = g_out_list[v1].first;
         i2 < g_out_list[v1].first + g_out_list[v1].second; ++i2) {
      const RawEdge &e2 = g_raw_edge[i2];
      const uint &v2 = e2.v;
      const ulong &w2 = e2.w;
      if (v2 <= st || !cmpWeight(w1, w2)) continue;
      for (uint i3 = g_out_list[v2].first;
           i3 < g_out_list[v2].first + g_out_list[v2].second; ++i3) {
        const RawEdge &e3 = g_raw_edge[i3];
        const uint &v3 = e3.v;
        const ulong &w3 = e3.w;
        if (v3 < st || v3 == v1 || !cmpWeight(w2, w3)) continue;
        if (v3 == st) {
          if (cmpWeight(w3, w1)) {
            // g_result[0][st].emplace_back(st, v1, v2);
            // g_answer_len[0][st] += (g_node_ids[st].len + g_node_ids[v1].len +
            // g_node_ids[v2].len);
            g_result[0][st].emplace_back(i1, i3);
            g_answer_len[0][st] += (e1.len + e3.half_len);
          }
          continue;
        }

        for (uint i4 = g_out_list[v3].first;
             i4 < g_out_list[v3].first + g_out_list[v3].second; ++i4) {
          const RawEdge &e4 = g_raw_edge[i4];
          const uint &v4 = e4.v;
          const ulong &w4 = e4.w;
          if (!(data.dist[v4] & 1)) continue;
          if (!cmpWeight(w3, w4))
            continue;
          else if (v4 == st) {
            if (cmpWeight(w4, w1)) {
              // g_result[1][st].emplace_back(st, v1, v2, v3);
              // g_answer_len[1][st] += (g_node_ids[st].len + g_node_ids[v1].len
              // + g_node_ids[v2].len + g_node_ids[v3].len);
              g_result[1][st].emplace_back(i1, i3);
              g_answer_len[1][st] += (e1.len + e3.len);
            }
            continue;
          } else if (v4 == v2 || v4 == v1)
            continue;
          for (uint i5 = g_out_list[v4].first;
               i5 < g_out_list[v4].first + g_out_list[v4].second; ++i5) {
            const RawEdge &e5 = g_raw_edge[i5];
            const uint &v5 = e5.v;
            const ulong &w5 = e5.w;
            if (!(data.dist[v5] & 2)) continue;
            if (!cmpWeight(w4, w5))
              continue;
            else if (v5 == st) {
              if (cmpWeight(w5, w1)) {  // TODO 这里好像不用判断了
                // g_result[2][st].emplace_back(st, v1, v2, v3, v4);
                // g_answer_len[2][st] += (g_node_ids[st].len +
                // g_node_ids[v1].len + g_node_ids[v2].len + g_node_ids[v3].len +
                // g_node_ids[v4].len);
                g_result[2][st].emplace_back(i1, i3, i5);
                g_answer_len[2][st] += (e1.len + e3.len + e5.half_len);
              }
              continue;
            } else if (v5 == v1 || v5 == v2 || v5 == v3)
              continue;
            for (uint i6 = g_out_list[v5].first;
                 i6 < g_out_list[v5].first + g_out_list[v5].second; ++i6) {
              const RawEdge &e6 = g_raw_edge[i6];
              const uint &v6 = e6.v;
              const ulong &w6 = e6.w;
              if (!(data.dist[v6] & 4)) continue;
              if (!cmpWeight(w5, w6))
                continue;
              else if (v6 == st) {
                if (cmpWeight(w6, w1)) {
                  // g_result[3][st].emplace_back(st, v1, v2, v3, v4, v5);
                  // g_answer_len[3][st] += (g_node_ids[st].len +
                  // g_node_ids[v1].len + g_node_ids[v2].len + g_node_ids[v3].len
                  // + g_node_ids[v4].len + g_node_ids[v5].len);
                  g_result[2][st].emplace_back(i1, i3, i5);
                  g_answer_len[2][st] += (e1.len + e3.len + e5.len);
                }
                continue;
              } else if (v6 == v1 || v6 == v2 || v6 == v3 || v6 == v4)
                continue;
              const ulong &w7 = data.end_weight[v6];
              if (cmpWeight(w6, w7) && cmpWeight(w7, w1)) {
                // g_result[4][st].emplace_back(st, v1, v2, v3, v4, v5, v6);
                // g_answer_len[4][st] += (g_node_ids[st].len +
                // g_node_ids[v1].len + g_node_ids[v2].len + g_node_ids[v3].len +
                // g_node_ids[v4].len + g_node_ids[v5].len + g_node_ids[v6].len);
                g_result[2][st].emplace_back(i1, i3, i5, i6);
                g_answer_len[2][st] +=
                    (e1.len + e3.len + e5.len + e6.len - e6.half_len);
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
}

void findCycleSkip(uint pid) {
  for (uint i = pid; i < g_node_num; i += kThreadNum) {
    if (g_in_list[i].empty() || !g_out_list[i].second) continue;
    // printf("find%d\n",i);
    findCycleAt(i, g_worker_data[pid]);
  }
}

void writePrePare() {
  int R = 0;
  for (auto &t : g_result) {
    for (auto &i : t) {
      R += i.size();
    }
  }
  cout << R << "cycles\n";

  // int *x;
  /*    char *temp_c = new char[20];

      ulong all_answer_len = 0;

      for (uint i = 0; i < 5; ++i) {
          for (uint j = 0; j < g_node_num; ++j) {
              all_answer_len += g_answer_len[i][j];
          }
      }

      sprintf(temp_c, "%d\n", R);
      all_answer_len += strlen(temp_c);


      int fd = open(g_predict_file.c_str(), O_RDWR | O_CREAT, 0666);

      fallocate(fd, 0, 0, all_answer_len);
      //ftruncate(fd,all_answer_len);

      char *answer_mmap =
              (char *) mmap(NULL, all_answer_len, PROT_READ | PROT_WRITE,
     MAP_SHARED, fd, 0);

      lseek(fd, all_answer_len - 1, SEEK_SET);
      write(fd, "\0", 1);
      close(fd);
      char *p = answer_mmap;
      for (int i = 0; i < strlen(temp_c); ++i) {
          *p = temp_c[i];
          p++;
      }
      g_path_mmap_begin = p;*/
}

void writeResult(int pid) {
  char *p = g_path_mmap_begin;
  for (uint i = 0; i < 5; ++i) {
    for (uint j = 0; j < g_node_num; ++j) {
      if (j % kThreadNum == pid) {
        for (auto &x : g_result[i][j]) {
          for (uint y = 0; y < i + 2; ++y) {
            // memcpy_tiny(p, g_node_ids[x.path[y]].val,
            // g_node_ids[x.path[y]].len); p += g_node_ids[x.path[y]].len;
          }
          // memcpy_tiny(p, g_node_ids[x.path[i + 2]].val, g_node_ids[x.path[i +
          // 2]].len);//TODO 这里可以-1 p += g_node_ids[x.path[i + 2]].len;
          *(p - 1) = '\n';
        }
      } else {
        p += g_answer_len[i][j];
      }
    }
  }
}

void findCycleTheading(uint pid) {
#ifdef LOCAL
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif

  //    g_worker_data[pid].result[0].reserve(kMaxResult / 512);
  //    g_worker_data[pid].result[1].reserve(kMaxResult / 128);
  //    g_worker_data[pid].result[2].reserve(kMaxResult / 32);
  //    g_worker_data[pid].result[3].reserve(kMaxResult / 8);
  //    g_worker_data[pid].result[4].reserve(kMaxResult / 2);

  // findCycleInRange(g_start_idx_list[pid], g_end_idx_list[pid], pid);
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
  exit(0);
  writeResult(pid);
#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("process %d write data time %fs\n", pid, t4 - t3);
#endif
}

void solve() {
  loadTestData();

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