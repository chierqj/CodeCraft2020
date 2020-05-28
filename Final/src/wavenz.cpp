#include <fcntl.h>
#include <pthread.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ext/pb_ds/priority_queue.hpp>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <vector>

#ifdef DEBUG
#include "timer.h"
#endif

using namespace std;

const char* test_data_file = "/data/test_data.txt";
const char* result_file = "/projects/student/result.txt";
// const char* result_file = "answer.txt";

const int nthread = 8;  // 线程数

const int maxn = 2500000;  // 最大节点数

bool dense = false;
// 数据读取buffer
uint32_t loadbuf[nthread][maxn][3];
uint32_t loadnum[nthread];

// 图数据结构（前向星）
struct Edge {
  uint32_t id;
  uint32_t weight;
};

Edge From[maxn];  // 正图
Edge To[maxn];    // 反图

uint32_t phead[maxn];  // 各id起始位置
uint32_t plen[maxn];   // 各id转账数
uint32_t qhead[maxn];
uint32_t qlen[maxn];

// 简易的固定长度vector，效率和数组相近
class Vector32 {
 public:
  Vector32() { _data = new uint32_t[maxn]; }
  void push_back(uint32_t t) { _data[_size++] = t; }
  int size() { return _size; }
  void clear() { _size = 0; }
  uint32_t* begin() { return _data; }
  uint32_t* end() { return _data + _size; }
  uint32_t& operator[](int index) { return _data[index]; }

 private:
  uint32_t* _data;
  int _size;
};
Vector32 ID;  //　真实ID表

// 简易的固定桶数hash_table，采用开放寻址法解决冲突
class Hash_map {
 public:
  Hash_map() { map.resize(m); }
  void insert(uint32_t k) {
    // 如果不存在则插入，否则直接返回，作用等同：
    // if(!Map.count(key)) Map.insert(key);
    int i = 0;
    int hash_val = 0;
    while (i < m) {
      hash_val = hash(k, i);
      if (map[hash_val].val == -1) {
        map[hash_val].key = k;
        map[hash_val].val = ID.size();
        map[hash_val].len = 1;
        ID.push_back(k);
        hashes.push_back(hash_val);
        cnt++;
        break;
      }
      if (map[hash_val].key == k) {
        map[hash_val].len++;
        break;
      } else
        i++;
    }
  }
  int search(uint32_t k) {
    // 搜索
    int i = 0, hash_val = 0;
    while (i < m) {
      hash_val = hash(k, i);
      if (map[hash_val].val == -1) break;
      if (map[hash_val].key == k)
        return map[hash_val].val;
      else
        i++;
    }
    return -1;
  }
  void sort_hash() {
    // 将hash值排序，确保id映射前后相对大小不变化
    sort(hashes.begin(), hashes.end(),
         [&](int a, int b) { return map[a].key < map[b].key; });
    for (int i = 0; i < hashes.size(); ++i) {
      map[hashes[i]].val = i;
    }
  }
  int size() { return cnt; }

 private:
  // 常数取质数，且 m1 < m
  const int m = 2222281;
  const int m1 = 2205167;
  int cnt = 0;
  // 数据结构
  struct data {
    uint32_t len = 0;
    uint32_t key;
    int val = -1;
  };
  vector<data> map;
  vector<int> hashes;
  uint32_t hash(uint32_t k, int i) {
    // 哈希函数
    return k % m + i;  // 一次哈希
    // return (k % m + i * (m1 - k % m1)) % m; // 双重哈希
  }
} Map;

struct Node {
  uint32_t dis;
  uint32_t id;
};
bool operator<(const Node& a, const Node& b) { return a.dis < b.dis; }
bool operator>(const Node& a, const Node& b) { return a.dis > b.dis; }
class prio_queue_parse {
 public:
  prio_queue_parse() {
    A = new Node[maxn];
    size = 0;
  }
  bool empty() { return size == 0; }
  void clear() { size = 0; }
  void push(Node node) {
    A[++size].dis = 0xffffffff;
    increase(size, node);
  }
  Node top() { return A[1]; }
  void pop() {
    A[1] = A[size];
    size--;
    heapify(1);
  }

 private:
  int parent(int curr) { return (curr >> 1); }
  int left(int curr) { return (curr << 1); }
  int right(int curr) { return ((curr << 1) + 1); }
  void heapify(int curr) {
    A[0] = A[curr];
    int l, r, largest;
    while (1) {
      l = left(curr), r = right(curr);
      if (l <= size && A[0] > A[l])
        largest = l;
      else
        largest = 0;
      if (r <= size && A[largest] > A[r]) largest = r;
      if (largest != 0) {
        A[curr] = A[largest];
        curr = largest;
      } else
        break;
    }
    A[curr] = A[0];
  }
  void increase(int curr, Node& target) {
    while (curr > 1 && A[parent(curr)] > target) {
      A[curr] = A[parent(curr)];
      curr = parent(curr);
    }
    A[curr] = target;
  }
  Node* A;
  size_t size;
};

class prio_queue_dense {
 public:
  prio_queue_dense() {
    A = new Node[maxn];
    index = new uint32_t[maxn]();
    size = 0;
  }
  bool empty() { return size == 0; }
  size_t getsize() { return size; }
  void clear() { size = 0; }
  void push(Node node) {
    int pos = index[node.id];
    if (pos == 0) {
      A[++size] = node;
      up(size);
    } else {
      A[pos].dis = node.dis;
      up(pos);
    }
  }
  Node top() { return A[1]; }
  void pop() {
    index[A[1].id] = 0;
    A[1] = A[size];
    index[A[size--].id] = 0;
    if (size) down(1);
  }

 private:
  int parent(int curr) { return (curr >> 1); }
  int left(int curr) { return (curr << 1); }
  int right(int curr) { return ((curr << 1) + 1); }
  void down(int curr) {
    A[0] = A[curr];
    int l, r, min;
    while (1) {
      l = left(curr), r = right(curr);
      if (l <= size && A[l] < A[0])
        min = l;
      else
        min = 0;
      if (r <= size && A[r] < A[min]) min = r;
      if (min != 0) {
        A[curr] = A[min];
        index[A[curr].id] = curr;
        curr = min;
      } else
        break;
    }
    A[curr] = A[0];
    index[A[curr].id] = curr;
  }
  void up(int curr) {
    Node temp = A[curr];
    while (curr > 1 && temp < A[parent(curr)]) {
      A[curr] = A[parent(curr)];
      index[A[curr].id] = curr;
      curr = parent(curr);
    }
    A[curr] = temp;
    index[temp.id] = curr;
  }
  Node* A;
  size_t size;
  uint32_t* index;  // 可以换成u16和u32各一个
};

class Prev {
 public:
  void push_back(const uint32_t& v, const uint32_t& u) {
    A[qhead[v] + len[v]++] = u;
  }
  void push(const uint32_t& v, const uint32_t& u) {
    len[v] = 0;
    A[qhead[v] + len[v]] = u;
    len[v] = 1;
  }
  uint32_t size(const uint32_t& v) { return len[v]; }
  void clear(const uint32_t& v) { len[v] = 0; }
  void clear_all() { memset(len, 0, ID.size() * sizeof(uint32_t)); }
  uint32_t* operator[](const uint32_t& v) { return A + qhead[v]; }

 private:
  uint32_t A[maxn];
  uint16_t len[maxn];  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
};

Prev back;
uint16_t In[maxn];
double res[nthread][maxn];
double restemp[nthread][maxn];
bool vis[nthread][maxn];
uint32_t dis[nthread][maxn];
uint16_t num[nthread][maxn];
Prev pre[nthread];
bool skip[maxn];

void go_back(uint32_t curr, const double& cent, int& tid) {
  res[tid][curr] += In[curr] * cent;
  for (int i = 0; i < back.size(curr); ++i) {
    go_back(back[curr][i], cent + 1, tid);
  }
}

atomic<int> cnt;
atomic<int> curNode;  // 原子计数，用于dfs负载均衡
void* find_thread_parse(void* arg) {
#ifdef DEBUG
  cout << __func__ << endl;
  Timer t;
#endif
  int tid = *(int*)arg;
  memset(dis[tid], 0xfe, ID.size() * sizeof(uint32_t));
  memset(num[tid], 0, ID.size() * sizeof(uint16_t));
  Vector32 clear;
  prio_queue_parse Q;
  while (1) {
    uint32_t s = curNode++;
    if (s >= ID.size()) break;
    if (skip[s]) continue;
#ifdef DEBUG
    if (s % 1000 == 0) cout << s << "/" << ID.size() << "   ";
    if (s % 1000 == 0) t.end();
#endif
    Q.clear();
    clear.clear();
    for (int i = phead[s]; i < phead[s + 1]; ++i) {
      num[tid][i] = 1;
    }
    num[tid][s] = 1;
    dis[tid][s] = 0;
    Q.push({dis[tid][s], s});
    while (!Q.empty()) {
      uint32_t u = Q.top().id;
      Q.pop();
      if (vis[tid][u]) continue;
      vis[tid][u] = true;
      clear.push_back(u);
      for (int i = phead[u]; i < phead[u + 1]; ++i) {
        uint32_t v = To[i].id;      // *
        if (vis[tid][v]) continue;  // *
        if (dis[tid][u] + To[i].weight < dis[tid][v]) {
          dis[tid][v] = dis[tid][u] + To[i].weight;
          num[tid][v] = num[tid][u];
          pre[tid].push(v, u);
          Q.push({dis[tid][v], v});
        } else if (dis[tid][u] + To[i].weight == dis[tid][v]) {
          num[tid][v] += num[tid][u];
          pre[tid].push_back(v, u);
        }
      }
    }
    for (int i = clear.size() - 1; i > 0; --i) {
      uint32_t& w = clear[i];
      double temp = (1 + restemp[tid][w]) / num[tid][w];
      for (int j = 0; j < pre[tid].size(w); ++j) {
        uint32_t& v = pre[tid][w][j];
        restemp[tid][v] += temp * num[tid][v];
      }
    }

    for (int i = 0; i < back.size(s); ++i) {
      go_back(back[s][i], restemp[tid][s] + 1, tid);
    }

    res[tid][s] += In[s] * restemp[tid][s];
    restemp[tid][s] = 0;
    vis[tid][s] = false;
    dis[tid][s] = 0xfefefefe;

    double temp = (1 + In[s]);
    for (int i = 1; i < clear.size(); ++i) {
      res[tid][clear[i]] += temp * restemp[tid][clear[i]];
      restemp[tid][clear[i]] = 0;
      vis[tid][clear[i]] = false;
      dis[tid][clear[i]] = 0xfefefefe;
    }
  }
}

void* find_thread_dense(void* arg) {
#ifdef DEBUG
  cout << __func__ << endl;
  Timer t;
#endif
  int tid = *(int*)arg;
  memset(dis[tid], 0xfe, ID.size() * sizeof(uint32_t));
  memset(num[tid], 0, ID.size() * sizeof(uint16_t));
  Vector32 clear;
  prio_queue_dense Q;
  while (1) {
    uint32_t s = curNode++;
    if (s >= ID.size()) break;
    if (skip[s]) continue;
    cnt++;
#ifdef DEBUG
    if (s % 1000 == 0) cout << s << "/" << ID.size() << "   ";
    if (s % 1000 == 0) t.end();
#endif
    Q.clear();
    clear.clear();
    for (int i = phead[s]; i < phead[s + 1]; ++i) {
      num[tid][i] = 1;
    }
    num[tid][s] = 1;
    dis[tid][s] = 0;
    Q.push({dis[tid][s], s});
    while (!Q.empty()) {
      uint32_t u = Q.top().id;
      Q.pop();
      if (vis[tid][u]) continue;
      vis[tid][u] = true;
      clear.push_back(u);
      for (int i = phead[u]; i < phead[u + 1]; ++i) {
        uint32_t v = To[i].id;      // *
        if (vis[tid][v]) continue;  // *
        if (dis[tid][u] + To[i].weight < dis[tid][v]) {
          dis[tid][v] = dis[tid][u] + To[i].weight;
          num[tid][v] = num[tid][u];
          Q.push({dis[tid][v], v});
        } else if (dis[tid][u] + To[i].weight == dis[tid][v]) {
          num[tid][v] += num[tid][u];
        }
      }
    }
    for (int i = clear.size() - 1; i > 0; --i) {
      uint32_t& w = clear[i];
      double temp = (1 + restemp[tid][w]) / num[tid][w];
      for (int j = qhead[w]; j < qhead[w + 1]; ++j) {
        uint32_t& v = From[j].id;
        if (dis[tid][v] + From[j].weight == dis[tid][w]) {
          restemp[tid][v] += temp * num[tid][v];
        }
      }
    }

    for (int i = 0; i < back.size(s); ++i) {
      go_back(back[s][i], restemp[tid][s] + 1, tid);
    }

    res[tid][s] += In[s] * restemp[tid][s];
    restemp[tid][s] = 0;
    vis[tid][s] = false;
    dis[tid][s] = 0xfefefefe;

    double temp = (1 + In[s]);
    for (int i = 1; i < clear.size(); ++i) {
      res[tid][clear[i]] += temp * restemp[tid][clear[i]];
      restemp[tid][clear[i]] = 0;
      vis[tid][clear[i]] = false;
      dis[tid][clear[i]] = 0xfefefefe;
    }
  }
}

void topo_sort() {
  queue<uint32_t> Q;
  for (int i = 0; i < ID.size(); ++i) {
    if (qlen[i] == 0) Q.push(i);
  }
  while (!Q.empty()) {
    uint32_t curr = Q.front();
    Q.pop();
    if (plen[curr] == 1) skip[curr] = true;
    for (int i = phead[curr]; i < phead[curr + 1]; ++i) {
      if (plen[curr] == 1) {
        In[To[i].id] += (In[curr] + 1);
        back.push_back(To[i].id, curr);
      }
      if (--qlen[To[i].id] == 0) Q.push(To[i].id);
    }
  }
}

void find_shortest() {
#ifdef DEBUG
  cout << __func__ << endl;
  Timer t;
#endif

  topo_sort();
  int nthread = 8;
  pthread_t threads[nthread];
  int tid[nthread];
  for (int i = 0; i < nthread; ++i) {
    tid[i] = i;
    if (dense)
      pthread_create(&threads[i], NULL, find_thread_dense, (void*)&tid[i]);
    else
      pthread_create(&threads[i], NULL, find_thread_parse, (void*)&tid[i]);
  }
  for (int i = 0; i < nthread; ++i) pthread_join(threads[i], NULL);
#ifdef DEBUG

#endif
}

void* sort_thread(void* arg) {
  // 对每个节点的邻接表排序：正向正序（保证结果为字典序），反向反序（用于反向dfs提前break退出）
  int tid = *(int*)arg;
  for (int i = tid; i < ID.size(); i += nthread) {
    sort(To + phead[i], To + phead[i] + plen[i],
         [](const Edge& a, const Edge& b) { return a.id < b.id; });
    sort(From + qhead[i], From + qhead[i] + qlen[i],
         [](const Edge& a, const Edge& b) { return a.id > b.id; });
  }
}

void* build_thread(void* arg) {
  // 构造前向星
  // 1.统计每个节点邻接点个数
  // 2.计算每个节点开始位置
  // 3.遍历所有节点构造前向星
  int tid = *(int*)arg;
  uint32_t from, to, weight;
  if (tid == 0) {  // 正图
    int* curlen = new int[ID.size() + 1]();
    for (int k = 0; k < nthread; ++k) {
      for (int i = 0; i < loadnum[k]; ++i) {
        from = loadbuf[k][i][0];
        plen[from]++;
      }
    }
    phead[0] = 0;
    for (int i = 1; i <= ID.size(); ++i) {
      phead[i] = phead[i - 1] + plen[i - 1];
    }
    for (int k = 0; k < nthread; ++k) {
      for (int i = 0; i < loadnum[k]; ++i) {
        from = loadbuf[k][i][0];
        to = loadbuf[k][i][1];
        To[phead[from] + curlen[from]].id = to;
        To[phead[from] + curlen[from]++].weight = loadbuf[k][i][2];
      }
    }
  } else {  // 反图
    int* curlen = new int[ID.size() + 1]();
    for (int k = 0; k < nthread; ++k) {
      for (int i = 0; i < loadnum[k]; ++i) {
        to = loadbuf[k][i][1];
        qlen[to]++;
      }
    }
    qhead[0] = 0;
    for (int i = 1; i <= ID.size(); ++i) {
      qhead[i] = qhead[i - 1] + qlen[i - 1];
    }
    for (int k = 0; k < nthread; ++k) {
      for (int i = 0; i < loadnum[k]; ++i) {
        to = loadbuf[k][i][1];
        from = loadbuf[k][i][0];
        From[qhead[to] + curlen[to]].id = from;
        From[qhead[to] + curlen[to]++].weight = loadbuf[k][i][2];
      }
    }
  }
}

char* file;
int file_size;
void* load_thread(void* args) {
  // 多线程读图
  int tid = *(int*)args;
  int size = file_size / nthread;
  if (tid == nthread - 1) size = file_size - (nthread - 1) * size;

  char* start = file + tid * (file_size / nthread);
  char* curr = start;

  // 确保两个线程不读到分割位置的同一行
  if (tid != 0 && *(curr - 1) != '\n')
    while (*curr++ != '\n')
      ;

  uint32_t from, to, weight, temp = 0;
  uint8_t state = 0;
  char ch;
  while (1) {
    ch = *curr;
    if (ch == ',') {
      state ? to = temp : from = temp;
      state = 1 - state;
      temp = 0;
    } else if (ch == '\r' || ch == '\n') {
      if (temp) {
        loadbuf[tid][loadnum[tid]][0] = from;
        loadbuf[tid][loadnum[tid]][1] = to;
        loadbuf[tid][loadnum[tid]][2] = temp;
        loadnum[tid]++;
      }

      if (ch == '\r') curr++;
      if (curr - start + 1 >= size) break;
      state = temp = 0;
    } else {
      temp = temp * 10 + ch - '0';
    }
    curr++;
  }
}
void* clear_thread(void* args) {
  // 删除记录中出度为０的节点
  // 并将记录id映射为{ 0 １ ２ ... n-1 }
  int tid = *(int*)args;
  uint32_t from, to;
  for (int i = 0; i < loadnum[tid]; ++i) {
    from = Map.search(loadbuf[tid][i][0]);
    to = Map.search(loadbuf[tid][i][1]);
    loadbuf[tid][i][0] = from;
    loadbuf[tid][i][1] = to;
  }
}

void* shit_thread(void* args) {
  // 哈希排序，确保映射前后相对大小不变
  int tid = *(int*)args;
  uint32_t from, to;
  if (tid == 0) {
    Map.sort_hash();
  } else {
    sort(ID.begin(), ID.end());
  }
}

void load_data() {
/*
    @brief: load data from file "/data/test_data.txt"
    @method: mmap
*/
#ifdef DEBUG
  cout << __func__ << endl;
  Timer t;
#endif
  // mmap
  struct stat statue;
  int fd = open(test_data_file, O_RDONLY);
  fstat(fd, &statue);
  file_size = statue.st_size;
  file = (char*)mmap(NULL, statue.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);
  // 1.多线程读数据到loadbuf中
  pthread_t threads[nthread];
  int tid[nthread];
  for (int i = 0; i < nthread; ++i) {
    tid[i] = i;
    pthread_create(&threads[i], NULL, load_thread, (void*)&tid[i]);
  }
  for (int i = 0; i < nthread; ++i) pthread_join(threads[i], NULL);

  // 2.哈希（仅映射 {u, v, w} 中的 u，因为不在 U 中出现的 id　肯定不成环）
  for (int k = 0; k < nthread; ++k) {
    for (int i = 0; i < loadnum[k]; ++i) {
      Map.insert(loadbuf[k][i][0]);
      Map.insert(loadbuf[k][i][1]);
    }
  }
  // 3.哈希排序，确保映射前后相对大小不变
  // for(int i = 0; i < 2; ++i){
  //     tid[i] = i;
  //     pthread_create(&threads[i], NULL, shit_thread, (void*)&tid[i]);
  // }
  // for(int i = 0; i < 2; ++i)
  //     pthread_join(threads[i], NULL);
  // 4.将不用的记录删除(使得后续不用再查询哈系表)
  for (int i = 0; i < nthread; ++i) {
    tid[i] = i;
    pthread_create(&threads[i], NULL, clear_thread, (void*)&tid[i]);
  }
  for (int i = 0; i < nthread; ++i) pthread_join(threads[i], NULL);

  // 5.构造前向星
  for (int i = 0; i < 2; ++i) {
    tid[i] = i;
    pthread_create(&threads[i], NULL, build_thread, (void*)&tid[i]);
  }
  for (int i = 0; i < 2; ++i) pthread_join(threads[i], NULL);

  // ６.前向星排序
  // for(int i = 0; i < nthread; ++i){
  //     tid[i] = i;
  //     pthread_create(&threads[i], NULL, sort_thread, (void*)&tid[i]);
  // }
  // for(int i = 0; i < nthread; ++i)
  //     pthread_join(threads[i], NULL);
}

int main(int argc, char** argv) {
#ifdef DEBUG
  Timer t;
#endif
  load_data();

  int num = 0;
  for (int i = 0; i < nthread; ++i) {
    num += loadnum[i];
  }
  if (num * 1.0 / ID.size() > 10) dense = true;

  find_shortest();
  for (int i = 0; i < ID.size(); ++i) {
    for (int j = 1; j < nthread; ++j) {
      res[0][i] += res[j][i];
    }
  }
  // for(int i = 0; i < ID.size(); ++i) cout << ID[i] << " ";
  // cout << endl;

  vector<int> index;
  for (int i = 0; i < ID.size(); ++i) {
    index.push_back(i);
  }
  sort(index.begin(), index.end(), [](int a, int b) {
    if (abs(res[0][a] - res[0][b]) <= 0.0001) return ID[a] < ID[b];
    return res[0][a] > res[0][b];
  });

  FILE* fp = fopen(result_file, "w");
  // sort(res, res + ID.size(), greater<double>());
  int sz = 100;
  if (ID.size() < 100) sz = ID.size();
  for (int i = 0; i < sz; ++i) {
    fprintf(fp, "%d,%.3lf\n", ID[index[i]], res[0][index[i]]);
  }
  fclose(fp);
#ifdef DEBUG
  cout << __func__ << endl;
#endif
  return 0;
}
