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
#include <sys/types.h>

#include <random>

using namespace std;
bool IS_LOCAL = false;  // run at local or submmit to server
/**
 * 本端账号ID和对端账号ID为一个32位的正整数
 * 转账金额为一个32位的正整数
 * 转账记录最多为300万条
 *
 * 账号A给账号B最多转账一次
 * **/

// TODO 指针换数组访问
// TODO 边找边写
// TODO 存0.2*w 3*w
// TODO 连逗号一起存 finish
// TODO vector存id  同时计算长度  最后一次性拷贝到mmap里
// TODO reserve
// TODO 不使用path了，直接把点存进去 然后判断%cyclelen==0输出\n
// TODO 存环的时候计算answerlen，然后分配空间mmap，直接拷贝到mmap分配的空间上
// TODO  多进程的节点分配为0 8 16 24   || 1 9 17 25 ||...
// TODO node相关的换成600w的数组

// typedef pair<uint32_t, uint32_t> Edge;
typedef struct Edge {
  uint32_t first;
  long long second;

  Edge(uint32_t x, long long w) : first(x), second(w){};
} Edge;
struct MyKeyHashHasher {
  inline size_t operator()(const int &k) const noexcept {
    return k & 0x7fffff;  // std::hash<int>{}(k.key);
  }
};
struct Path {
  int path[7];  //不初始化 初始化浪费时间

  Path() = default;

  Path(int a, int b, int c) {
    path[0] = a;
    path[1] = b;
    path[2] = c;
  }

  Path(int a, int b, int c, int d) {
    path[0] = a;
    path[1] = b;
    path[2] = c;
    path[3] = d;
  }

  Path(int a, int b, int c, int d, int e) {
    path[0] = a;
    path[1] = b;
    path[2] = c;
    path[3] = d;
    path[4] = e;
  }

  Path(int a, int b, int c, int d, int e, int f) {
    path[0] = a;
    path[1] = b;
    path[2] = c;
    path[3] = d;
    path[4] = e;
    path[5] = f;
  }

  Path(int a, int b, int c, int d, int e, int f, int g) {
    path[0] = a;
    path[1] = b;
    path[2] = c;
    path[3] = d;
    path[4] = e;
    path[5] = f;
    path[6] = g;
  }
};

const int MAX_RESULT = 20000000;
const int MAX_N = 6000000;
const int T = 8;                  //进程/线程数量
const int MAX_IDVAL_LENGTH = 11;  //整数表示的ID有多长
int RAND_KEY;                     //随机值用于shmget
int RAND_VISITED_FLG;

int N = 0;  // node num
string test_file;
string predict_file;

vector<uint32_t> node_id;
// vector<int> node_weight;
vector<char *> node_ids;       //节点id的char*,int compare快,char* write快
vector<size_t> node_ids_size;  //记录char*的大小
vector<vector<Edge>> in_list;
vector<vector<Edge>> out_list;
unordered_map<uint32_t, uint32_t> id_index_map;
vector<uint32_t> arg_list;  // arg_sort结果
uint32_t worker_answer_len[5];
vector<Path> worker_result[5];
// int worker_result_num[5];

// id升序
bool cmpAdj(const Edge &a, const Edge &b);

// id降序
bool cmpAdjReverse(const Edge &a, const Edge &b);

// char*转int，将从*p到*q的一段char*内存转换成int/longlong
inline void myScan(char *p, const char *q, uint32_t &x);
inline void myScan(char *p, const char *q, long long &x);
// x对y转账w元
inline void addEdge(char *&x_str, char *&y_str, uint32_t x, uint32_t y,
                    size_t x_len, size_t y_len, uint32_t w);

//解析输入文件
void analyzeBuffMmap(char *buffer, int MAXS);

//读取输入文件并建立邻接表
bool loadTestData();

//节点id argsort
vector<uint32_t> idArgSort();

//找从一个起点开始的环
void findCycleAt(uint32_t st, int *dist, long long *end_weight);

//找一个区间内所有点的环
void findCycleInRange(uint32_t st, uint32_t ed);

//进程/线程找环
void findCycleProcess(int pid);

//进程同步
void waitOtherProcess(int pid, int *shared_buffer, int target_mode);

//写结果
void writeResultProcess(int pid);

void solve();

/***********************************************************************************************/

bool cmpAdj(const Edge &a, const Edge &b) {
  return node_id[a.first] < node_id[b.first];
}

bool cmpAdjReverse(const Edge &a, const Edge &b) {
  return node_id[a.first] > node_id[b.first];
}

inline void myScan(char *p, const char *q, uint32_t &x) {
  x = *p++ - '0';
  while (p != q) {
    x = (x << 3) + (x << 1) + (*p - '0');
    ++p;
  }
}
inline void myScan(char *p, const char *q, long long &x) {
  x = *p++ - '0';
  while (p != q) {
    x = (x << 3) + (x << 1) + (*p - '0');
    ++p;
  }
}

inline void addEdge(char *&x_str, char *&y_str, uint32_t x, uint32_t y,
                    size_t x_len, size_t y_len, uint32_t w) {
  int xid, yid;
  // TODO 用默认值==0替代//存.find结果减少寻址
  auto it_x = id_index_map.find(x);
  auto it_y = id_index_map.find(y);
  auto it_ed = id_index_map.end();

  if (it_x == it_ed) {
    node_id.emplace_back(x);
    node_ids.emplace_back(x_str);
    node_ids_size.emplace_back(x_len);
    in_list.emplace_back();
    out_list.emplace_back();
    id_index_map[x] = N;
    xid = N;
    N++;
    it_ed = id_index_map.end();
  } else {
    xid = it_x->second;
  }

  if (it_y == it_ed) {
    node_id.emplace_back(y);
    node_ids.emplace_back(y_str);
    node_ids_size.emplace_back(y_len);
    in_list.emplace_back();
    out_list.emplace_back();
    id_index_map[y] = N;
    yid = N;
    N++;
  } else {
    yid = it_y->second;
  }

  out_list[xid].emplace_back(Edge(yid, w));
  in_list[yid].emplace_back(Edge(xid, w));
}

void analyzeBuffMmap(char *buffer, int MAXS) {
  //预分配空间
  id_index_map.reserve(MAX_N);
  node_id.reserve(MAX_N);
  node_ids.reserve(MAX_N);
  node_ids_size.reserve(MAX_N);
  in_list.reserve(MAX_N);
  out_list.reserve(MAX_N);

  const char split = ',';
  const char split_line = '\r';  // TODO 线上数据为/r/n
  char *p = buffer;
  char *q = buffer;
  // string xs, ys;
  uint32_t x, y;
  long long w;
  char *xs, *ys;
  size_t lx, ly;  // length

  while (q - buffer < MAXS) {
    while (*q != split) q++;
    xs = p;
    lx = q - p + 1;

    myScan(p, q, x);
    p = ++q;
    while (*q != split) q++;

    ys = p;
    ly = q - p + 1;
    myScan(p, q, y);
    p = ++q;

    while (*q != split_line) q++;
    myScan(p, q, w);
    q += 2;
    p = q;

    addEdge(xs, ys, x, y, lx, ly, w);
  }
  for (int i = 0; i < N; ++i) {
    if (in_list[i].size() > 1)
      sort(in_list[i].begin(), in_list[i].end(), cmpAdjReverse);
    if (out_list[i].size() > 1)
      sort(out_list[i].begin(), out_list[i].end(), cmpAdj);
  }
}

bool loadTestData() {
  ifstream fin(test_file.c_str(), std::ios::binary);
  if (!fin) {
    cout << "打开训练文件失败" << endl;
    exit(0);
  }
  if (IS_LOCAL) {
    clock_t start, end;
    start = clock();
    int MAX_BUFFER_SIZE = fin.seekg(0, std::ios::end).tellg();
    char *file_buffer = new char[MAX_BUFFER_SIZE];  // vector<char>(MAXS);
    fin.seekg(0, std::ios::beg).read(file_buffer, MAX_BUFFER_SIZE);
    analyzeBuffMmap(file_buffer, MAX_BUFFER_SIZE);
    end = clock();
    cout << N << "nodes\n";
    cout << "Load data time : " << ((double)end - start) / CLOCKS_PER_SEC
         << "s\n";
  } else {
    int MAX_BUFFER_SIZE = fin.seekg(0, std::ios::end).tellg();
    char *file_buffer = new char[MAX_BUFFER_SIZE];  // vector<char>(MAXS);
    fin.seekg(0, std::ios::beg).read(file_buffer, MAX_BUFFER_SIZE);
    analyzeBuffMmap(file_buffer, MAX_BUFFER_SIZE);
  }
  fin.close();
  return true;
}

vector<uint32_t> idArgSort() {
  vector<uint32_t> array_index(N, 0);
  for (uint32_t i = 0; i < N; ++i) array_index[i] = i;

  std::sort(array_index.begin(), array_index.end(),
            [](int pos1, int pos2) { return (node_id[pos1] < node_id[pos2]); });

  return array_index;
}

void findCycleAt(uint32_t st, int *dist, long long *end_weight) {
  uint32_t v1, v2, v3, v4, v5, v6, i1, i2, i3, i4, i5, i6;
  long long w1, w2, w3, w4, w5, w6, w7;
  uint32_t st_id = node_id[st];
  vector<int> changed_dist;
  // changed_dist.reserve(2500);
  // TODO w*3 w*5可能溢出  用longlong或者预先存的时候 判断>INTMAX/3 INTMAX/5
  // 若大于则直接置为INTMAX
  for (i1 = 0; i1 < in_list[st].size(); ++i1) {
    v1 = in_list[st][i1].first;
    if (node_id[v1] < st_id) break;
    w1 = in_list[st][i1].second;
    dist[v1] = 1;
    end_weight[v1] = w1;
    changed_dist.push_back(v1);
    for (i2 = 0; i2 < in_list[v1].size(); ++i2) {
      v2 = in_list[v1][i2].first;
      if (node_id[v2] <= st_id) break;
      w2 = in_list[v1][i2].second;  // X=w2  Y=w1  X<=5Y  3X>=Y
      if (w2 > 5 * w1 || 3 * w2 < w1) continue;
      dist[v2] = dist[v2] < 2 ? 1 : 2;
      changed_dist.push_back(v2);
      for (i3 = 0; i3 < in_list[v2].size(); ++i3) {
        v3 = in_list[v2][i3].first;
        if (node_id[v3] <= st_id) break;
        // else if (c == a)
        // continue;
        w3 = in_list[v2][i3].second;
        if (w3 > 5 * w2 || 3 * w3 < w2) continue;
        dist[v3] = dist[v3] < 3 ? dist[v3] : 3;
        changed_dist.push_back(v3);
      }
    }
  }
  dist[st] = 0;
  Path trace{};
  trace.path[0] = st;
  // TODO 提前计算size st~v2
  for (i1 = 0; i1 < out_list[st].size(); ++i1) {
    v1 = out_list[st][i1].first;
    if (node_id[v1] < st_id) continue;
    w1 = out_list[st][i1].second;
    trace.path[1] = v1;
    for (i2 = 0; i2 < out_list[v1].size(); ++i2) {
      v2 = out_list[v1][i2].first;
      if (node_id[v2] < st_id || v2 == st) continue;
      w2 = out_list[v1][i2].second;
      if (w1 > 5 * w2 || 3 * w1 < w2) continue;
      trace.path[2] = v2;
      for (i3 = 0; i3 < out_list[v2].size(); ++i3) {
        v3 = out_list[v2][i3].first;
        if (node_id[v3] < st_id) continue;
        w3 = out_list[v2][i3].second;
        if (w2 > 5 * w3 || 3 * w2 < w3) continue;
        if (v3 == v1) continue;
        if (v3 == st) {
          // addAnswer(st, v1, v2);
          if (w3 <= 5 * w1 && 3 * w3 >= w1) {
            worker_result[0].emplace_back(trace);
            worker_answer_len[0] +=
                (node_ids_size[st] + node_ids_size[v1] + node_ids_size[v2]);
          }
          continue;
        }

        // worker_prefix_size = 0;
        trace.path[3] = v3;
        for (i4 = 0; i4 < out_list[v3].size(); ++i4) {
          v4 = out_list[v3][i4].first;
          if (dist[v4] > 3 || node_id[v4] < st_id) continue;
          w4 = out_list[v3][i4].second;
          if (w3 > 5 * w4 || 3 * w3 < w4) continue;
          if (v4 == v2 || v4 == v1) continue;
          if (v4 == st) {
            if (w4 <= 5 * w1 && 3 * w4 >= w1) {
              worker_result[1].emplace_back(trace);
              worker_answer_len[1] += (node_ids_size[st] + node_ids_size[v1] +
                                       node_ids_size[v2] + node_ids_size[v3]);
            }
            continue;
          }
          trace.path[4] = v4;
          for (i5 = 0; i5 < out_list[v4].size(); ++i5) {
            v5 = out_list[v4][i5].first;
            if (dist[v5] > 2 || node_id[v5] < st_id) continue;
            w5 = out_list[v4][i5].second;
            if (w4 > 5 * w5 || 3 * w4 < w5) continue;
            if (v5 == v1 || v5 == v2 || v5 == v3) continue;
            if (v5 == st) {
              if (w5 <= 5 * w1 && 3 * w5 >= w1) {  // TODO 这里好像不用判断了
                worker_result[2].emplace_back(trace);
                worker_answer_len[2] +=
                    (node_ids_size[st] + node_ids_size[v1] + node_ids_size[v2] +
                     node_ids_size[v3] + node_ids_size[v4]);
              }
              continue;
            }
            trace.path[5] = v5;
            for (i6 = 0; i6 < out_list[v5].size(); ++i6) {
              v6 = out_list[v5][i6].first;
              if (dist[v6] > 1 || node_id[v6] < st_id) {
                continue;
              }
              w6 = out_list[v5][i6].second;
              if (w5 > 5 * w6 || 3 * w5 < w6) continue;
              if (v6 == v1 || v6 == v2 || v6 == v3 || v6 == v4) continue;

              if (v6 == st) {
                if (w6 <= 5 * w1 && 3 * w6 >= w1) {
                  worker_result[3].emplace_back(trace);
                  worker_answer_len[3] +=
                      (node_ids_size[st] + node_ids_size[v1] +
                       node_ids_size[v2] + node_ids_size[v3] +
                       node_ids_size[v4] + node_ids_size[v5]);
                }
                continue;
              }

              trace.path[6] = v6;
              // addAnswer(worker_path_prefix, worker_prefix_size, st, v1, v2,
              // v3, v4, v5, v6);
              w7 = end_weight[v6];
              if (w7 <= 5 * w1 && 3 * w7 >= w1) {
                worker_result[4].emplace_back(trace);
                worker_answer_len[4] +=
                    (node_ids_size[st] + node_ids_size[v1] + node_ids_size[v2] +
                     node_ids_size[v3] + node_ids_size[v4] + node_ids_size[v5] +
                     node_ids_size[v6]);
              }
            }
          }
        }
      }
    }
  }
  for (auto x : changed_dist) {
    dist[x] = 4;
  }
  dist[st] = 4;
}

void findCycleInRange(uint32_t st, uint32_t ed) {
  int *dist = new int[N];
  long long *adj_weight = new long long[N];
  for (int i = 0; i < N; ++i) {
    // inTrace[i] = false;
    dist[i] = 4;
  }
  for (unsigned int &i : worker_answer_len) {
    i = 0;
  }
  for (int i = st; i < ed; ++i) {
    if (in_list[arg_list[i]].empty() || out_list[arg_list[i]].empty()) continue;
    findCycleAt(arg_list[i], dist, adj_weight);
  }
}

void findCycleProcess(int pid) {
  // TODO 自适应大小
  int step = N / 54;
  vector<int> st_list = {0,         2 * step,  4 * step,  7 * step,
                         11 * step, 17 * step, 25 * step, 35 * step};
  vector<int> ed_list = {2 * step,  4 * step,  7 * step,  11 * step,
                         17 * step, 25 * step, 35 * step, N};
  findCycleInRange(st_list[pid], ed_list[pid]);
}

void waitOtherProcess(int pid, int *shared_buffer, int target_mode) {
  shared_buffer[pid] = target_mode;
  bool flg;
  while (true) {
    flg = true;
    for (int i = 0; i < T; ++i) {
      if (shared_buffer[i] != target_mode) flg = false;
    }
    if (flg) {
      break;
    } else
      usleep(10);
  }
}

void writeResultProcess(int pid) {
  int shmid = shmget(RAND_KEY, 4096, IPC_CREAT | 0666);
  int *shared_buffer = (int *)shmat(shmid, NULL, 0);
  int WRITE_STATE = T * 4;
  /* 0~T用于状态同步  RAND_VISITED_FLG表示找完环了
   * RAND_VISITED_FLG-1表示写完数据了 T~2T表示各自的环的数量
   * 4T用于记录轮到谁写了  mmap应该用不到
   * 5T~9T记录各自的answer length*/

  int R = 0;
  for (int j = 0; j < 5; ++j) {
    R += worker_result[j].size();
    // worker_answer_len[j] = worker_answer_p[j] - worker_answer[j];
  }

  shared_buffer[WRITE_STATE] = 0;
  shared_buffer[pid + T] = R;
  for (int i = 0; i < 5; ++i) {
    shared_buffer[(i + 5) * T + pid] = worker_answer_len[i];
  }

  //开始同步
  waitOtherProcess(pid, shared_buffer, RAND_VISITED_FLG);

  int all_R = 0;
  int all_answer_len = 0;
  for (int i = 0; i < T; ++i) {
    all_R += shared_buffer[i + T];
  }

  char *temp_c = new char[20];
  sprintf(temp_c, "%d\n", all_R);
  all_answer_len += strlen(temp_c);
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < T; ++j) {
      all_answer_len += shared_buffer[(i + 5) * T + j];
    }
  }

  int fd = open(predict_file.c_str(), O_RDWR | O_CREAT, 0666);

  char *answer_mmap = (char *)mmap(NULL, all_answer_len, PROT_READ | PROT_WRITE,
                                   MAP_SHARED, fd, 0);

  //主进程写总数
  if (pid == T - 1) {
    fallocate(fd, 0, 0, all_answer_len);  // TODO 放在pid==0里
    printf("%d cycles\n", all_R);
    // printf("all answer len %d\n",all_answer_len);
    lseek(fd, all_answer_len - 1, SEEK_SET);
    write(fd, "\0", 1);

    char *p = answer_mmap;
    for (int i = 0; i < strlen(temp_c); ++i) {
      *p = temp_c[i];
      p++;
    }
  }
  close(fd);

  //写答案
  char *p = answer_mmap + strlen(temp_c);
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < T; ++j) {
      usleep(1);  // TODO 不加这个会卡死 为什么？
      if (j == pid) {
        for (auto &x : worker_result[i]) {
          for (int y = 0; y < i + 2; ++y) {
            memcpy(p, node_ids[x.path[y]], node_ids_size[x.path[y]]);
            p += node_ids_size[x.path[y]];
          }
          memcpy(p, node_ids[x.path[i + 2]],
                 node_ids_size[x.path[i + 2]]);  // TODO 这里可以-1
          p += node_ids_size[x.path[i + 2]];
          *(p - 1) = '\n';
        }
      } else {
        p += shared_buffer[(i + 5) * T + j];
      }
    }
  }

  //结束同步
  waitOtherProcess(pid, shared_buffer, RAND_VISITED_FLG + 1);
  // printf("process %d end\n",pid);
  exit(0);
}

void solve() {
  loadTestData();
  default_random_engine e(time(0));
  uniform_int_distribution<int> u(0, 1008611);
  RAND_KEY = u(e);
  RAND_VISITED_FLG = u(e);

  arg_list = idArgSort();
  pid_t Process[T] = {0};
  int pid = 0;
  for (int i = 1; i < T; ++i) {
    if (pid == 0) {
      int x = fork();
      Process[i] = x;
      if (x == 0) {
        pid = i;
        break;
      }
    }
  }
  // printf("this is process %d\n",pid);
  findCycleProcess(pid);
  writeResultProcess(pid);
}

int main(int argc, char *argv[]) {
  if (argc == 2 && strcmp(argv[1], "local") == 0) IS_LOCAL = true;

  if (!IS_LOCAL) {
    test_file = "./data/big/test_data.txt";
    predict_file = "./data/big/result.txt";
    // solve();
    exit(0);
  } else {
    cout << "local test" << endl;

    vector<string> all_testFile = {"testdata/standard.txt"};
    // all_testFile={"testdata/1004812/test_data.txt"};
    vector<string> all_predictFile = {"data/result.txt"};

    for (int i = 0; i < all_testFile.size(); i += 1) {
      cout << "Example " << i << endl;
      test_file = all_testFile[i];
      predict_file = all_predictFile[i];
      solve();
    }
  }

  return 0;
}
