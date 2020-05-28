#include <bits/stdc++.h>
#include <fcntl.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
using namespace std;

struct timeval t1, t2;
double timeuse;

typedef pair<uint64_t, uint32_t> P1;
typedef pair<uint32_t, uint32_t> P2;
uint32_t ALLDATASIZE = 0;
uint32_t MAXDATASIZE = 1000000;
uint32_t THX = 8;
uint64_t MAXz = 0;
atomic<uint32_t> RW(1);

uint32_t outys[1000000];
unordered_map<uint32_t, uint32_t> in_ys;

uint32_t DIL[7500000];
uint32_t DATA[5000000];

uint32_t DLEN[1000000];
typedef struct {
  bool vis[1000000];
  uint16_t pre[1000000];
  uint8_t lpe[1000000];
  uint32_t mpe[1000000][50];
  uint32_t sta[1000000];
  double tpp[1000000];

} THD;

THD thdx[8];

uint64_t dis64x[8][1000000];
uint32_t dis32x[8][1000000];
double ans[8][1000000];

string testfile = "/data/test_data.txt";
string out_file = "/projects/student/result.txt";
// string testfile="data/data5/test_data.txt";
// string out_file="data/my_out1.txt";
struct node64 {
  uint64_t first;
  uint32_t second;
  node64(const uint64_t &_first, const uint32_t &_second)
      : first(_first), second(_second) {}
  bool operator<(const node64 &r) const { return first > r.first; }
};
struct node32 {
  uint32_t first;
  uint32_t second;
  node32() {
    first = 0;
    second = 0;
  };
  node32(const uint32_t &_first, const uint32_t &_second)
      : first(_first), second(_second) {}
  bool operator<(const node32 &r) const { return first < r.first; }
  node32 operator=(const node32 &r) {
    first = r.first;
    second = r.second;
    return *this;
  }
};
class HEAP {
 public:
  node32 node[2000000];
  // int id[1000000];
  int size;
  bool empty() { return size == 0 ? 1 : 0; }
  inline void swapn(node32 &x, node32 &y) {
    auto p = x;
    x = y;
    y = p;
  }

  void up(node32 x) {
    // node[id[x.second]].first=x.first;
    for (int i = size, j = i >> 1; j; i = j, j >>= 1) {
      if (node[i] < node[j]) {
        swapn(node[i], node[j]);
        // swap(id[node[i].second],id[node[j].second]);
      } else
        break;
    }
    return;
  }
  void push(node32 x) {
    node[++size] = x;
    // id[x.second]=size;
    up(x);
  }
  node32 top() { return node[1]; }

  void pop() {
    // id[node[size].second]=1;
    // id[1]=0;
    swapn(node[1], node[size]);
    size--;
    int i = 1;
    for (int j = 2; j <= size; i = j, j <<= 1) {
      if (j < size && node[j + 1] < node[j]) j++;
      if (node[j] < node[i]) {
        swapn(node[i], node[j]);
        // swap(id[node[i].second],id[node[j].second]);
      } else
        break;
    }
  }
};

HEAP heapx[8];
void th_in_ys(int index, int lenend) {
  for (int i = index; i < lenend; i += 3) {
    DIL[i] = in_ys.at(DIL[i]);
    DIL[i + 1] = in_ys.at(DIL[i + 1]);
  }
}

void getdata() {
  int fd = open(testfile.c_str(), O_RDONLY);
  uint64_t len = lseek(fd, 0, SEEK_END);
  char *buf = (char *)mmap(NULL, len, PROT_READ, MAP_PRIVATE, fd, 0);

  vector<uint32_t> paixu;

  char *p = buf;
  uint32_t tmp1 = 0;
  uint32_t tmp2 = 0;
  uint32_t key = 0;
  uint32_t Y = 0;
  char *km = buf;
  km += len;
  while (p < km) {
    while (*p != ',') {
      tmp1 = int(*p - '0') + tmp1 * 10;
      p++;
    }
    p++;
    while (*p != ',') {
      tmp2 = int(*p - '0') + tmp2 * 10;
      p++;
    }
    p++;
    while (*p != '\n' && *p != '\r') {
      key = int(*p - '0') + key * 10;
      p++;
    }
    if (*p == '\r') p++;
    if (key != 0) {
      paixu.push_back(tmp1);
      paixu.push_back(tmp2);
      DIL[Y++] = tmp1;
      DIL[Y++] = tmp2;
      DIL[Y++] = key;
      MAXz += key;
    }
    tmp1 = 0;
    tmp2 = 0;
    key = 0;
    p++;
  }

  munmap(buf, len);
  close(fd);

  sort(paixu.begin(), paixu.end());  // haoshi

  ALLDATASIZE = unique(paixu.begin(), paixu.end()) - paixu.begin();

  for (int i = 1; i <= ALLDATASIZE; i++) {
    in_ys.insert({paixu[i - 1], i});
    outys[i] = paixu[i - 1];
  }
  int canshu[9] = {0};
  for (int i = 1; i < 8; i++) canshu[i] = (int)(Y / 24) * 3 * i;

  canshu[8] = Y;
  vector<thread *> thr;
  for (int i = 0; i < 8; i++)
    thr.push_back(new thread(th_in_ys, canshu[i], canshu[i + 1]));
  for (auto th : thr) th->join();

  vector<vector<P2>> data_g(ALLDATASIZE + 1);

  for (int i = 0; i < Y; i += 3)

    data_g[DIL[i]].push_back(make_pair(DIL[i + 1], DIL[i + 2]));

  for (auto &line : data_g) sort(line.begin(), line.end());
  int num = 0;
  DLEN[0] = 0;
  for (int i = 0; i <= ALLDATASIZE; i++) {
    for (auto x : data_g[i]) {
      DATA[num++] = x.first;
      DATA[num++] = x.second;
      // DKEY[num++]=x.second;
    }
    DLEN[i + 1] = num;
  }
}

bool cmpsize(pair<double, uint32_t> &line1, pair<double, uint32_t> &line2) {
  if (abs(line1.first - line2.first) > 0.0001) return line1.first > line2.first;
  return line1.second < line2.second;
}

void out_data() {
  for (int i = 1; i < THX; i++) {
    for (int j = 1; j <= ALLDATASIZE; j++) ans[0][j] += ans[i][j];
  }

  vector<pair<double, uint32_t>> out_data(ALLDATASIZE);
  for (int i = 1; i <= ALLDATASIZE; i++)
    out_data[i - 1] = pair<double, uint32_t>(ans[0][i], i);
  sort(out_data.begin(), out_data.end(), cmpsize);
  int MINXX = ALLDATASIZE > 100 ? 100 : ALLDATASIZE;
  FILE *fp = fopen(out_file.c_str(), "w");
  for (int i = 0; i < MINXX; i++) {
    fprintf(fp, "%u,%.3f\n", outys[out_data[i].second], out_data[i].first);
  }
}
void dij(uint32_t s, uint32_t index, THD &thdx, uint64_t *dis64,
         uint64_t &INF) {
  priority_queue<node64> que;  // greater less
  int nM = 0;
  dis64[s] = 0;
  thdx.pre[s] = 1;
  que.push(node64(0, s));
  while (!que.empty()) {
    node64 p = que.top();
    que.pop();
    auto &v = p.second;
    if (thdx.vis[v]) continue;
    thdx.vis[v] = 1;
    thdx.sta[nM++] = v;
    for (uint32_t i = DLEN[v]; i < DLEN[v + 1]; i += 2)  //遍历与顶点v相连的边
    {
      if (thdx.vis[DATA[i]] == 0) {
        auto &w = DATA[i];
        auto KK = dis64[v] + DATA[i + 1];
        if (dis64[w] > KK) {
          thdx.lpe[w] = 0;
          dis64[w] = KK;
          thdx.pre[w] = thdx.pre[v];
          thdx.mpe[w][thdx.lpe[w]++] = v;
          que.push(node64(KK, w));
        } else if (dis64[w] == KK) {
          thdx.mpe[w][thdx.lpe[w]++] = v;
          thdx.pre[w] += thdx.pre[v];
        }
      }
    }
  }

  dis64[s] = INF;
  thdx.tpp[s] = 0;
  thdx.lpe[s] = 0;
  thdx.vis[s] = 0;
  while (nM) {
    nM--;
    auto &w = thdx.sta[nM];

    for (auto i = 0; i < thdx.lpe[w]; i++) {
      auto &v = thdx.mpe[w][i];
      if (v != s)
        thdx.tpp[v] += (1.0 + thdx.tpp[w]) * thdx.pre[v] / thdx.pre[w];
    }

    ans[index][w] += thdx.tpp[w];
    dis64[w] = INF;
    thdx.tpp[w] = 0;
    thdx.vis[w] = 0;
    thdx.lpe[w] = 0;
    thdx.pre[w] = 0;
  }
}

void dij1(uint32_t s, uint32_t index, THD &thdx, uint32_t *dis32, uint32_t &INF,
          HEAP &heap) {
  // priority_queue<node32> que;//greater less

  int nM = 0;
  dis32[s] = 0;
  thdx.pre[s] = 1;
  heap.push(node32(0, s));
  while (heap.size) {
    node32 p = heap.node[1];
    heap.pop();
    auto &v = p.second;
    if (thdx.vis[v]) {
      continue;
    }
    thdx.vis[v] = 1;
    thdx.sta[nM++] = v;
    for (uint32_t i = DLEN[v]; i < DLEN[v + 1]; i += 2)  //遍历与顶点v相连的边
    {
      if (thdx.vis[DATA[i]] == 0) {
        auto &w = DATA[i];
        auto KK = dis32[v] + DATA[i + 1];
        if (dis32[w] > KK) {
          thdx.lpe[w] = 0;
          dis32[w] = KK;
          thdx.pre[w] = thdx.pre[v];
          thdx.mpe[w][thdx.lpe[w]++] = v;
          // if(heap.id[w])
          //   heap.up(node32(KK,w));
          // else
          heap.push(node32(KK, w));

        } else if (dis32[w] == KK) {
          thdx.mpe[w][thdx.lpe[w]++] = v;
          thdx.pre[w] += thdx.pre[v];
        }
      }
    }
  }

  dis32[s] = INF;
  thdx.tpp[s] = 0;
  thdx.lpe[s] = 0;
  thdx.vis[s] = 0;
  while (nM) {
    nM--;
    auto &w = thdx.sta[nM];

    for (auto i = 0; i < thdx.lpe[w]; i++) {
      auto &v = thdx.mpe[w][i];
      if (v != s)
        thdx.tpp[v] += (1.0 + thdx.tpp[w]) * thdx.pre[v] / thdx.pre[w];
    }

    ans[index][w] += thdx.tpp[w];
    dis32[w] = INF;
    thdx.tpp[w] = 0;
    thdx.vis[w] = 0;
    thdx.lpe[w] = 0;
    thdx.pre[w] = 0;
  }
}

void thrs1(uint32_t index) {
  uint64_t INF64 = pow(2, 55);
  uint32_t INF32 = pow(2, 32);
  if (MAXz > INF32) {
    THD &thd = thdx[index];
    for (int j = 0; j <= ALLDATASIZE; j++) {
      thd.vis[j] = 0;
      thd.pre[j] = 0;
      thd.lpe[j] = 0;
      thd.tpp[j] = 0;
      ans[index][j] = 0;
      dis64x[index][j] = INF64;
    }
    uint32_t i = RW.fetch_add(1);
    while (i <= ALLDATASIZE) {
      dij(i, index, thd, dis64x[index], INF64);
      if (i % 10000 == 0) {
        cout << i << "/" << ALLDATASIZE << " ";
        gettimeofday(&t2, NULL);
        timeuse = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec) / 1000000.0;
        printf(" %d dij Use Time:%fs\n", i, timeuse);
      }
      i = RW.fetch_add(1);
    }

  } else {
    THD &thd = thdx[index];
    HEAP &heap = heapx[index];
    heap.size = 0;
    for (int j = 0; j <= ALLDATASIZE; j++) {
      thd.vis[j] = 0;
      thd.pre[j] = 0;
      thd.lpe[j] = 0;
      thd.tpp[j] = 0;
      // heap.id[j]=0;

      ans[index][j] = 0;
      dis32x[index][j] = INF32;
    }
    uint32_t i = RW.fetch_add(1);
    while (i <= ALLDATASIZE) {
      dij1(i, index, thd, dis32x[index], INF32, heap);
      if (i % 10000 == 0) {
        cout << i << "/" << ALLDATASIZE << " ";
        gettimeofday(&t2, NULL);
        timeuse = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec) / 1000000.0;
        printf(" %d dij Use Time:%fs\n", i, timeuse);
      }
      i = RW.fetch_add(1);
    }
  }
}

int main() {
  gettimeofday(&t1, NULL);

  getdata();

  vector<thread *> thr;

  for (int k = 0; k < THX; k++) thr.push_back(new thread(thrs1, k));

  for (auto &th : thr) th->join();

  out_data();

  gettimeofday(&t2, NULL);
  timeuse = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec) / 1000000.0;
  printf(" end Use Time:%fs\n", timeuse);

  exit(0);

  return 0;
}