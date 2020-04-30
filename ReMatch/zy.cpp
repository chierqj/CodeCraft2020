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

typedef pair<uint32_t, uint32_t> Edge;
const int MAX_RESULT = 20000000;

const int MAX_IDVAL_LENGTH = 11;  //整数表示的ID有多长

int N = 0;  // node num
string test_file;
string predict_file;

vector<uint32_t> node_id;
vector<int> node_weight;
vector<char *> node_ids;
vector<size_t> node_ids_size;
vector<vector<Edge>> in_list;
vector<vector<Edge>> out_list;

unordered_map<uint32_t, uint32_t> id_index_map;

// int node_id[MAX_N];//节点id
// const char *node_ids[MAX_N];//节点id的char*,int compare快,char* write快
// size_t node_ids_size[MAX_N];//记录char*的大小
// int in_list[MAX_N][MAX_DEG];//邻接表入边list
// int out_list[MAX_N][MAX_DEG];//邻接表出边list
// int in_list_size[MAX_N] = {0};//邻接表入边list大小
// int out_list_size[MAX_N] = {0};//邻接表出边list大小
// unordered_map<int, vector<pair<int, int>>> prePath[MAX_N];//前驱路径
// prePath[st][v]表示从v到st的所有长度3的合法路径,起点终点已知故只用存中间 int
// id_index_map[MAX_ID_RANGE];//线上数据范围如果过大要换成hashset

vector<uint32_t> arg_list;  // arg_sort结果
// int all_answer_len;

// vector<Path> all_result[5];
// const int T=4;
const int T = 8;

//对应T个线程的结果路径
// char *worker_answer[5];//结果路径保存为char*形式
// char *worker_answer_p[5];
uint32_t worker_answer_len[5];
// int worker_result_num[5];

struct Path {
  int path[7];

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

vector<Path> worker_result[5];

bool cmpAdj(const Edge &a, const Edge &b) {
  return node_id[a.first] < node_id[b.first];
}

bool cmpAdjReverse(const Edge &a, const Edge &b) {
  return node_id[a.first] > node_id[b.first];
}

// x对y转账
inline void addEdge(char *&xs, char *&ys, uint32_t x, uint32_t y, size_t lx,
                    size_t ly, uint32_t w) {
  int xid, yid;
  // TODO 用默认值==0替代//存.find结果减少寻址
  if (id_index_map.find(x) == id_index_map.end()) {
    node_id.push_back(x);
    node_ids.push_back(xs);
    node_ids_size.push_back(lx);
    in_list.emplace_back();
    out_list.emplace_back();
    id_index_map[x] = N;
    xid = N;
    N++;
  } else {
    xid = id_index_map[x];
  }

  if (id_index_map.find(y) == id_index_map.end()) {
    node_id.push_back(y);
    node_ids.push_back(ys);
    node_ids_size.push_back(ly);
    in_list.emplace_back();
    out_list.emplace_back();
    id_index_map[y] = N;
    yid = N;
    N++;
  } else {
    yid = id_index_map[y];
  }

  out_list[xid].push_back({yid, w});
  in_list[yid].push_back({xid, w});
}

//将从*p到*q的一段char*内存转换成int
inline void myScan(char *p, const char *q, uint32_t &x) {
  x = *p++ - '0';
  while (p != q) {
    x = (x << 3) + (x << 1) + (*p - '0');
    ++p;
  }
}

//解析输入文件
void analyzeBuffMmap(char *buffer, int MAXS) {
  const char split = ',';
  const char split_line = '\n';
  char *p = buffer;
  char *q = buffer;
  // string xs, ys;
  uint32_t x, y, w;
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
    p = ++q;
    addEdge(xs, ys, x, y, lx, ly, w);
  }
  for (int i = 0; i < N; ++i) {
    if (in_list[i].size() > 1)
      sort(in_list[i].begin(), in_list[i].end(), cmpAdjReverse);
    if (out_list[i].size() > 1)
      sort(out_list[i].begin(), out_list[i].end(), cmpAdj);
  }
}

//读取输入文件并建立邻接表
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

//对index按id排序
vector<uint32_t> idArgSort() {
  vector<uint32_t> array_index(N, 0);
  for (uint32_t i = 0; i < N; ++i) array_index[i] = i;

  std::sort(array_index.begin(), array_index.end(),
            [](int pos1, int pos2) { return (node_id[pos1] < node_id[pos2]); });

  return array_index;
}

void distFindCyclePathWorker(uint32_t st, int *dist, char *path_prefix) {
  uint32_t v1, v2, v3, v4, v5, v6, i1, i2, i3, i4, i5, i6;
  uint32_t st_id = node_id[st];
  vector<int> changed_dist;
  // changed_dist.reserve(2500);
  //    char *worker_path_prefix = path_prefix;
  //    uint32_t worker_prefix_size = 0;

  for (i1 = 0; i1 < in_list[st].size(); ++i1) {
    v1 = in_list[st][i1].first;
    if (node_id[v1] < st_id) break;
    dist[v1] = 1;
    changed_dist.push_back(v1);
    for (i2 = 0; i2 < in_list[v1].size(); ++i2) {
      v2 = in_list[v1][i2].first;
      if (node_id[v2] <= st_id) break;
      dist[v2] = dist[v2] < 2 ? 1 : 2;
      changed_dist.push_back(v2);
      for (i3 = 0; i3 < in_list[v2].size(); ++i3) {
        v3 = in_list[v2][i3].first;
        if (node_id[v3] <= st_id) break;
        //                else if (c == a)
        //                    continue;
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

    trace.path[1] = v1;
    for (i2 = 0; i2 < out_list[v1].size(); ++i2) {
      v2 = out_list[v1][i2].first;
      if (node_id[v2] < st_id || v2 == st) continue;
      trace.path[2] = v2;
      for (i3 = 0; i3 < out_list[v2].size(); ++i3) {
        v3 = out_list[v2][i3].first;
        if (node_id[v3] < st_id) continue;
        if (v3 == v1) continue;
        if (v3 == st) {
          // addAnswer(st, v1, v2);
          worker_result[0].emplace_back(trace);
          worker_answer_len[0] +=
              (node_ids_size[st] + node_ids_size[v1] + node_ids_size[v2]);
          continue;
        }

        // worker_prefix_size = 0;
        trace.path[3] = v3;
        for (i4 = 0; i4 < out_list[v3].size(); ++i4) {
          v4 = out_list[v3][i4].first;
          if (dist[v4] > 3 || node_id[v4] < st_id) continue;
          if (v4 == v2 || v4 == v1) continue;
          if (v4 == st) {
            // addAnswer(worker_path_prefix, worker_prefix_size, st, v1, v2,
            // v3);
            worker_result[1].emplace_back(trace);
            worker_answer_len[1] += (node_ids_size[st] + node_ids_size[v1] +
                                     node_ids_size[v2] + node_ids_size[v3]);
            continue;
          }
          trace.path[4] = v4;
          for (i5 = 0; i5 < out_list[v4].size(); ++i5) {
            v5 = out_list[v4][i5].first;
            if (dist[v5] > 2 || node_id[v5] < st_id) continue;
            if (v5 == v1 || v5 == v2 || v5 == v3) continue;
            if (v5 == st) {
              // addAnswer(worker_path_prefix, worker_prefix_size, st, v1, v2,
              // v3, v4);
              worker_result[2].emplace_back(trace);
              worker_answer_len[2] +=
                  (node_ids_size[st] + node_ids_size[v1] + node_ids_size[v2] +
                   node_ids_size[v3] + node_ids_size[v4]);
              continue;
            }
            trace.path[5] = v5;
            for (i6 = 0; i6 < out_list[v5].size(); ++i6) {
              v6 = out_list[v5][i6].first;
              if (dist[v6] > 1 || node_id[v6] < st_id) {
                continue;
              }
              if (v6 == st) {
                // addAnswer(worker_path_prefix, worker_prefix_size, st, v1, v2,
                // v3, v4, v5);
                worker_result[3].emplace_back(trace);
                worker_answer_len[3] +=
                    (node_ids_size[st] + node_ids_size[v1] + node_ids_size[v2] +
                     node_ids_size[v3] + node_ids_size[v4] + node_ids_size[v5]);
                continue;
              }
              if (v6 == v1 || v6 == v2 || v6 == v3 || v6 == v4) continue;
              trace.path[6] = v6;
              // addAnswer(worker_path_prefix, worker_prefix_size, st, v1, v2,
              // v3, v4, v5, v6);
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
  for (auto x : changed_dist) {
    dist[x] = 4;
  }
  dist[st] = 4;
}

void findCycleWorker(uint32_t st, uint32_t ed) {
  // printf("st%d ed%d\n",st,ed);
  // bool *inTrace = new bool[N];
  int *dist = new int[N];
  char *path_prefix = new char[MAX_IDVAL_LENGTH * 3];
  for (int i = 0; i < N; ++i) {
    // inTrace[i] = false;
    dist[i] = 4;
  }
  int decay_rate[5] = {100, 30, 10, 3, 1};
  for (int i = 0; i < 5; ++i) {
    //        worker_answer[i] = new char[MAX_RESULT * (i + 3) *
    //        MAX_IDVAL_LENGTH / decay_rate[i]]; worker_answer_p[i] =
    //        worker_answer[i];
    // worker_result_num[i] = 0;
    worker_answer_len[i] = 0;
  }
  for (int i = st; i < ed; ++i) {
    if (in_list[arg_list[i]].size() == 0 || out_list[arg_list[i]].size() == 0)
      continue;
    distFindCyclePathWorker(arg_list[i], dist, path_prefix);
  }

  //    if (pid==3)
  //        sleep(15);
}

void findCycleThreading(int pid) {
  int step = N / 200000;
  //    vector<int> st_list = {0, 4990 * step, 5001 * step, 12560 * step};
  //    vector<int> ed_list{4990 * step, 5001 * step, 12560 * step, N};
  vector<int> st_list = {0,
                         6250 * step,
                         9980 * step,
                         10001 * step,
                         11002 * step,
                         17000 * step,
                         25000 * step,
                         28000 * step};
  vector<int> ed_list{6250 * step,  9980 * step,  10001 * step, 11002 * step,
                      17000 * step, 25000 * step, 28000 * step, N};

  if (N < 200000 || N > 210000) {
    step = N / 54;
    st_list = {0,         2 * step,  4 * step,  7 * step,
               11 * step, 17 * step, 25 * step, 35 * step};
    ed_list = {2 * step,  4 * step,  7 * step,  11 * step,
               17 * step, 25 * step, 35 * step, N};
  }

  findCycleWorker(st_list[pid], ed_list[pid]);
}

int RAND_KEY;
int RAND_VISITED_FLG;

void waitWorkerResult(int pid) {
  int shmid = shmget(RAND_KEY, 4096, IPC_CREAT | 0666);
  int *addr = (int *)shmat(shmid, NULL, 0);
  int WRITE_STATE = T * 4;
  // 0~T用于状态同步  RAND_VISITED_FLG表示找完环了
  // RAND_VISITED_FLG-1表示写完数据了 T~2T表示各自的环的数量
  // 4T用于记录轮到谁写了 mmap应该用不到 5T~9T记录各自的answerlen

  int R = 0;

  for (int j = 0; j < 5; ++j) {
    R += worker_result[j].size();
    // worker_answer_len[j] = worker_answer_p[j] - worker_answer[j];
  }

  addr[WRITE_STATE] = 0;
  addr[pid + T] = R;
  for (int i = 0; i < 5; ++i) {
    addr[(i + 5) * T + pid] = worker_answer_len[i];
  }

  addr[pid] = RAND_VISITED_FLG;

  // printf("process %d begin\n",pid);
  //开始同步
  bool flg;
  while (true) {
    flg = true;
    for (int i = 0; i < T; ++i) {
      if (addr[i] != RAND_VISITED_FLG) flg = false;
    }
    if (flg) {
      break;
    } else
      usleep(10);
  }
  // printf("process %d write\n",pid);

  int all_R = 0;
  int all_answer_len = 0;
  for (int i = 0; i < T; ++i) {
    all_R += addr[i + T];
  }

  char *temp_c = new char[20];
  sprintf(temp_c, "%d\n", all_R);
  all_answer_len += strlen(temp_c);
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < T; ++j) {
      all_answer_len += addr[(i + 5) * T + j];
    }
  }

  int fd = open(predict_file.c_str(), O_RDWR | O_CREAT, 0666);
  // fallocate(fd, 0, 0, all_answer_len);
  char *answer_mmap = (char *)mmap(NULL, all_answer_len, PROT_READ | PROT_WRITE,
                                   MAP_SHARED, fd, 0);

  //主进程写总数据
  if (pid == 0) {
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
  // printf("process %d ready\n",pid);
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
        p += addr[(i + 5) * T + j];
      }

      // printf("%d %d %d\n", pid, i, j);
    }
  }

  // printf("process %d finish\n",pid);

  //结束同步
  addr[pid] = RAND_VISITED_FLG + 1;
  while (true) {
    flg = true;
    for (int i = 0; i < T; ++i) {
      if (addr[i] != RAND_VISITED_FLG + 1) flg = false;
    }
    if (flg) {
      break;
    } else
      usleep(10);
  }

  // printf("process %d end\n",pid);
  exit(0);
}

void solve() {
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
  findCycleThreading(pid);
  waitWorkerResult(pid);
}

int main(int argc, char *argv[]) {
  if (argc == 2 && strcmp(argv[1], "local") == 0) IS_LOCAL = true;

  if (!IS_LOCAL) {
    string testFile = "./data/1004812/test_data.txt";
    string predictFile = "./data/1004812/result.txt";
    test_file = testFile;
    predict_file = predictFile;
    loadTestData();
    solve();
    return 0;
  } else {
    cout << "local test" << endl;

    vector<string> all_testFile = {"./data/1004812/test_data.txt"};
    vector<string> all_predictFile = {"./data/1004812/result.txt"};

    for (int i = 0; i < all_testFile.size(); i += 1) {
      cout << "Example " << i << endl;
      test_file = all_testFile[i];
      predict_file = all_predictFile[i];
      loadTestData();
      solve();
    }
  }

  return 0;
}
