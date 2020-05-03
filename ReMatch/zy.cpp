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
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/time.h>
#include <sys/types.h>

#include <random>

using namespace std;
#define uint uint32_t
#define ulong uint64_t

string g_test_file;
string g_predict_file;
const uint kMaxResult = 20000008;
const uint kMaxE = 3000008;
const uint kMaxN = 6000008;
const uint kThreadNum = 4;     //进程/线程数量
const uint kMaxIdStrLen = 11;  //整数表示的ID有多长

typedef pair<uint, ulong> Edge;
/*struct Edge {
    ulong w;
    uint v;

    Edge() = default;

    Edge(uint x, ulong weight) : v(x), w(weight) {};
};*/

typedef struct RawEdge {
  ulong w;
  uint u;
  uint v;

  RawEdge() = default;

  RawEdge(uint x, uint y, ulong weight) : u(x), v(y), w(weight){};
} RawEdge;

typedef struct NodeIdstr {
  uint len;
  const char *val;

  NodeIdstr() = default;

  NodeIdstr(const char *p, uint l) : val(p), len(l){};
} Ids;

typedef struct CirclePath {
  int path[7];  //不初始化 初始化浪费时间

  CirclePath() = default;

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

  CirclePath(int a, int b, int c, int d, int e) {
    path[0] = a;
    path[1] = b;
    path[2] = c;
    path[3] = d;
    path[4] = e;
  }

  CirclePath(int a, int b, int c, int d, int e, int f) {
    path[0] = a;
    path[1] = b;
    path[2] = c;
    path[3] = d;
    path[4] = e;
    path[5] = f;
  }

  CirclePath(int a, int b, int c, int d, int e, int f, int g) {
    path[0] = a;
    path[1] = b;
    path[2] = c;
    path[3] = d;
    path[4] = e;
    path[5] = f;
    path[6] = g;
  }
} Path;

typedef struct WorkerData {
  uint changed_num = 0;
  uint answer_len[5]{0};
  ulong end_weight[kMaxN];
  uint changed_list[kMaxN];
  char dist[kMaxN];
  vector<Path> result[5];
  bool is_finish = false;
} WorkerData;

uint g_edge_num = 0;
uint g_node_num = 0;

vector<uint> g_start_idx_list;
vector<uint> g_end_idx_list;

uint g_node_id[kMaxN];
Ids g_unordered_node_ids[kMaxN];
Ids g_node_ids[kMaxN];

RawEdge g_raw_edge[kMaxE];
vector<Edge> g_in_list[kMaxN];
vector<Edge> g_out_list[kMaxN];

WorkerData g_worker_data[kThreadNum];

template <typename T>
inline void parseInteger(char *p, const char *q, T &x) {
  x = *p++ - '0';
  while (p != q) {
    x = (x << 3) + (x << 1) + (*p - '0');
    ++p;
  }
}

void addEdge(const char *&x_str, const char *&y_str, uint x_len, uint y_len,
             uint x, uint y, ulong w, unordered_map<uint, uint> &id_index_map) {
  // TODO 用默认值==0替代//存.find结果减少寻址
  auto it_x = id_index_map.find(x);
  auto it_y = id_index_map.find(y);
  auto it_ed = id_index_map.end();

  if (it_x == it_ed) {
    g_node_id[g_node_num] = x;
    g_unordered_node_ids[g_node_num].val = x_str;
    g_unordered_node_ids[g_node_num].len = x_len;

    // id_index_map[x] = N;
    id_index_map.insert(make_pair(x, g_node_num));
    ++g_node_num;
    it_ed = id_index_map.end();
  }

  if (it_y == it_ed) {
    g_node_id[g_node_num] = y;
    g_unordered_node_ids[g_node_num].val = y_str;
    g_unordered_node_ids[g_node_num].len = y_len;
    // id_index_map[y] = N;
    id_index_map.insert(make_pair(y, g_node_num));
    ++g_node_num;
  }
  g_raw_edge[g_edge_num].u = x;
  g_raw_edge[g_edge_num].v = y;
  g_raw_edge[g_edge_num].w = w;
  ++g_edge_num;
}

uint *idArgSort() {
  uint *arg_list = new uint[g_node_num];
  uint *arg_reflact = new uint[g_node_num];
  for (uint i = 0; i < g_node_num; ++i) arg_list[i] = i;

  std::sort(arg_list, arg_list + g_node_num, [](int pos1, int pos2) {
    return (g_node_id[pos1] < g_node_id[pos2]);
  });
  for (uint i = 0; i < g_node_num; ++i) {
    arg_reflact[arg_list[i]] = i;
    g_node_ids[i] = g_unordered_node_ids[arg_list[i]];
  }
  return arg_reflact;
}

void buildGraph(const uint *arg_reflect,
                unordered_map<uint, uint> &id_index_map, uint pid) {
  for (uint i = 0; i < g_edge_num; ++i) {
    RawEdge &e = g_raw_edge[i];
    if (e.u % kThreadNum == pid) {
      ulong w = e.w;
      uint x_id = arg_reflect[id_index_map[e.u]];
      uint y_id = arg_reflect[id_index_map[e.v]];
      g_out_list[x_id].emplace_back(Edge(y_id, w));
    }
    if (e.v % kThreadNum == pid) {
      ulong w = e.w;
      uint x_id = arg_reflect[id_index_map[e.u]];
      uint y_id = arg_reflect[id_index_map[e.v]];
      g_in_list[y_id].emplace_back(Edge(x_id, w));
    }
  }
}

bool cmpAdj(const Edge &a, const Edge &b) { return a.first < b.first; }

bool cmpAdjReverse(const Edge &a, const Edge &b) { return a.first > b.first; }

void sortAdjList(uint st, uint ed) {
  for (int i = st; i < ed; ++i) {
    if (g_in_list[i].size() > 1)
      sort(g_in_list[i].begin(), g_in_list[i].end(), cmpAdjReverse);
    if (g_out_list[i].size() > 1)
      sort(g_out_list[i].begin(), g_out_list[i].end(), cmpAdj);
  }
}

void analyzeBuff(char *buffer, uint max_buffer_size) {
  const char split = ',';
  const char split_line = '\n';  // TODO 线上数据为/r/n
  char *p = buffer;
  char *q = buffer;
  // string xs, ys;
  uint x, y, x_len, y_len;
  ulong w;
  const char *x_str, *y_str;
  unordered_map<uint, uint> id_index_map;
  id_index_map.reserve(1000000);
  id_index_map.rehash(1000000);

  while (q - buffer < max_buffer_size) {
    while (*q != split) q++;
    x_str = p;
    x_len = q - p + 1;

    parseInteger(p, q, x);
    p = ++q;
    while (*q != split) q++;

    y_str = p;
    y_len = q - p + 1;
    parseInteger(p, q, y);
    p = ++q;

    while (*q != '\r' && *q != '\n') q++;
    parseInteger(p, q, w);
    if (*q == '\r') ++q;
    p = ++q;
    addEdge(x_str, y_str, x_len, y_len, x, y, w, id_index_map);
  }
  cout << g_edge_num << " " << id_index_map.size() << endl;
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

  /*    for (int i=0;i<3;++i) {
          for (int k = 0; k < g_out_list[i].size(); ++k) {
              uint j=g_out_list[i][k].first;
              char a[20],b[20];
              memcpy(a,g_node_ids[i].val,g_node_ids[i].len-1);
              a[g_node_ids[i].len-1]='\0';
              memcpy(b,g_node_ids[j].val,g_node_ids[j].len-1);
              b[g_node_ids[j].len-1]='\0';
              printf("%d-%d : %s-%s %d \n",i,j,a,b,g_out_list[i][k].second);
          }
      }*/
}

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

inline double Pow4(double x) { return pow(x, 3.5); }

void setProcessNodeNum() {
  g_start_idx_list = {0};
  g_end_idx_list.clear();
  double ND = g_node_num;
  double total_cnt = 0;
  for (int i = 0; i < g_node_num; i += 25) {
    total_cnt += Pow4(g_out_list[i].size() * (1 - i / ND));
  }

  double step = total_cnt / kThreadNum;
  double next_bar = step;
  double cnt = 0;

  for (int i = 0; i < g_node_num; i += 25) {
    cnt += Pow4(g_out_list[i].size() * (1 - i / ND));
    ;
    if (cnt > next_bar) {
      g_start_idx_list.push_back(i + 1);
      g_end_idx_list.push_back(i + 1);
      next_bar += step;
    }
  }
  while (g_start_idx_list.size() < g_node_num) {
    g_start_idx_list.push_back(g_node_num);
    g_end_idx_list.push_back(g_node_num);
  }
  g_end_idx_list.push_back(g_node_num);
  for (int i = 0; i < kThreadNum; i++) printf("%d ", g_start_idx_list[i]);
  printf("%d\n", g_end_idx_list[kThreadNum - 1]);
}

inline bool cmpWeight(const ulong &x, const ulong &y) {
  return (x <= (y + (y << 2))) && ((x + (x << 1)) >= y);
}

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
  for (const auto &e1 : g_out_list[st]) {
    const uint &v1 = e1.first;
    const ulong &w1 = e1.second;
    if (v1 < st) continue;
    for (const auto &e2 : g_out_list[v1]) {
      const uint &v2 = e2.first;
      const ulong &w2 = e2.second;
      if (v2 <= st || !cmpWeight(w1, w2)) continue;
      for (const auto &e3 : g_out_list[v2]) {
        const uint &v3 = e3.first;
        const ulong &w3 = e3.second;
        if (v3 < st || v3 == v1 || !cmpWeight(w2, w3)) continue;
        if (v3 == st) {
          if (cmpWeight(w3, w1)) {
            data.result[0].emplace_back(st, v1, v2);
            data.answer_len[0] +=
                (g_node_ids[st].len + g_node_ids[v1].len + g_node_ids[v2].len);
          }
          continue;
        }

        for (const auto &e4 : g_out_list[v3]) {
          const uint &v4 = e4.first;
          if (!(data.dist[v4] & 1)) continue;
          if (!cmpWeight(w3, e4.second))
            continue;
          else if (v4 == st) {
            if (cmpWeight(e4.second, w1)) {
              data.result[1].emplace_back(st, v1, v2, v3);
              data.answer_len[1] += (g_node_ids[st].len + g_node_ids[v1].len +
                                     g_node_ids[v2].len + g_node_ids[v3].len);
            }
            continue;
          } else if (v4 == v2 || v4 == v1)
            continue;
          for (const auto &e5 : g_out_list[v4]) {
            const uint &v5 = e5.first;
            ;
            if (!(data.dist[v5] & 2)) continue;
            if (!cmpWeight(e4.second, e5.second))
              continue;
            else if (v5 == st) {
              if (cmpWeight(e5.second, w1)) {  // TODO 这里好像不用判断了
                data.result[2].emplace_back(st, v1, v2, v3, v4);
                data.answer_len[2] += (g_node_ids[st].len + g_node_ids[v1].len +
                                       g_node_ids[v2].len + g_node_ids[v3].len +
                                       g_node_ids[v4].len);
              }
              continue;
            } else if (v5 == v1 || v5 == v2 || v5 == v3)
              continue;
            for (const auto &e6 : g_out_list[v5]) {
              const uint &v6 = e6.first;
              if (!(data.dist[v6] & 4)) continue;
              if (!cmpWeight(e5.second, e6.second))
                continue;
              else if (v6 == st) {
                if (cmpWeight(e6.second, w1)) {
                  data.result[3].emplace_back(st, v1, v2, v3, v4, v5);
                  data.answer_len[3] +=
                      (g_node_ids[st].len + g_node_ids[v1].len +
                       g_node_ids[v2].len + g_node_ids[v3].len +
                       g_node_ids[v4].len + g_node_ids[v5].len);
                }
                continue;
              } else if (v6 == v1 || v6 == v2 || v6 == v3 || v6 == v4)
                continue;
              const ulong &w7 = data.end_weight[v6];
              if (cmpWeight(e6.second, w7) && cmpWeight(w7, w1)) {
                data.result[4].emplace_back(st, v1, v2, v3, v4, v5, v6);
                data.answer_len[4] += (g_node_ids[st].len + g_node_ids[v1].len +
                                       g_node_ids[v2].len + g_node_ids[v3].len +
                                       g_node_ids[v4].len + g_node_ids[v5].len +
                                       g_node_ids[v6].len);
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

void findCycleInRange(uint st, uint ed, uint pid) {
  st = st < g_node_num ? st : g_node_num;
  ed = ed < g_node_num ? ed : g_node_num;

  for (int i = st; i < ed; ++i) {
    if (g_in_list[i].empty() || g_out_list[i].empty()) continue;
    findCycleAt(i, g_worker_data[pid]);
  }
}

char *path_mmap_begin;
void writePrePare() {
  int R = 0;
  for (auto &t : g_worker_data) {
    for (auto &i : t.result) {
      R += i.size();
    }
  }
  cout << R << "cycles\n";

  // int *x;
  char *temp_c = new char[20];

  ulong all_answer_len = 0;

  for (auto &t : g_worker_data) {
    for (unsigned int i : t.answer_len) all_answer_len += i;
  }

  sprintf(temp_c, "%d\n", R);
  all_answer_len += strlen(temp_c);

  int fd = open(g_predict_file.c_str(), O_RDWR | O_CREAT, 0666);

  fallocate(fd, 0, 0, all_answer_len);

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
  path_mmap_begin = p;
}

void writeResult(int pid) {
  char *p = path_mmap_begin;
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < kThreadNum; ++j) {
      usleep(1);  // TODO 不加这个会卡死 为什么？
      if (j == pid) {
        for (auto &x : g_worker_data[j].result[i]) {
          for (int y = 0; y < i + 2; ++y) {
            memcpy(p, g_node_ids[x.path[y]].val, g_node_ids[x.path[y]].len);
            p += g_node_ids[x.path[y]].len;
          }
          memcpy(p, g_node_ids[x.path[i + 2]].val,
                 g_node_ids[x.path[i + 2]].len);  // TODO 这里可以-1
          p += g_node_ids[x.path[i + 2]].len;
          *(p - 1) = '\n';
        }
      } else {
        p += g_worker_data[j].answer_len[i];  // shared_buffer[(i + 5) * T + j];
      }
    }
  }
}

void findCycleTheading(uint pid) {
  // TODO 自适应大小
  //    uint32_t step = N / 54;
  //    start_node_id_list = {0, 2 * step, 4 * step, 7 * step, 11 * step, 17 *
  //    step, 25 * step, 35 * step}; end_node_id_list = {2 * step, 4 * step, 7 *
  //    step, 11 * step, 17 * step, 25 * step, 35 * step, N};
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  findCycleInRange(g_start_idx_list[pid], g_end_idx_list[pid], pid);
  gettimeofday(&tim, nullptr);
  double t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  cout << pid << " find cycle time " << t2 - t1 << "s\n";

  if (pid == 0) {
    while (true) {
      bool flg = true;
      for (uint i = 1; i < kThreadNum; i++) {
        flg = flg && g_worker_data[i].is_finish;
      }
      if (flg) break;
      usleep(10);
    }
    writePrePare();
  }
  g_worker_data[pid].is_finish = true;

  while (true) {
    bool flg = true;
    for (uint i = 0; i < kThreadNum; i++) {
      flg = flg && g_worker_data[i].is_finish;
    }
    if (flg) break;
    usleep(10);
  }
  writeResult(pid);
  gettimeofday(&tim, nullptr);
  double t3 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  cout << pid << " write data time " << t3 - t2 << "s\n";
}

void solve() {
  loadTestData();
  setProcessNodeNum();

  thread pool[kThreadNum - 1];
  for (uint i = 0; i < kThreadNum - 1; ++i) {
    pool[i] = thread(findCycleTheading, i);
  }
  findCycleTheading(kThreadNum - 1);
  for (auto &th : pool) th.join();
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