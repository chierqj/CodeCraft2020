#include <bits/stdc++.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
using namespace std;

//#define TEST

#define Inline __inline__ __attribute__((always_inline))

typedef uint8_t ui8;
typedef uint16_t ui16;
typedef uint32_t ui32;
typedef uint64_t ui64;

const ui8 THREAD_NUM = 8;
ui32 g_node_cnt;
char* g_buf;
struct Bi {
  ui32 tid;
  ui32 l;
  ui32 r;
};
vector<ui32> g_unhash_id;
vector<ui32> g_node_data[THREAD_NUM];
vector<ui32> g_edge_data[THREAD_NUM];

Inline ui32 MyAtoi(char* s) {
  ui32 ans = 0;
  while (47 < *s && *s < 58) {
    ans *= 10;
    ans += (*s - 48);
    s++;
  }
  return ans;
}

void ReadTestData(Bi& bi) {
  char* bg = g_buf + bi.l;
  const char* ed = g_buf + bi.r;
  char tmp[16];
  g_node_data[bi.tid].reserve(1000000);
  g_edge_data[bi.tid].reserve(500000);
  ui32 line[3];
  while (bg < ed) {
    for (ui32& ii : line) {
      ui8 index = 0;
      while (*(bg + index) != ',' && *(bg + index) != '\n') index++;
      memset(tmp, ' ', 16);
      memcpy(tmp, bg, 16);
      ii = MyAtoi(tmp);
      bg += (index + 1);
    }
    if (line[2] != 0) {
      g_node_data[bi.tid].emplace_back(line[0]);
      g_node_data[bi.tid].emplace_back(line[1]);
      g_edge_data[bi.tid].emplace_back(line[2]);
    }
  }
}

void Read(const string& file_name) {
#ifdef TEST
  auto start = std::chrono::steady_clock::now();
#endif
  int fd = open(file_name.c_str(), O_RDONLY);
  ui32 length = lseek(fd, 0, SEEK_END);
  g_buf = (char*)mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0);
  thread threads[THREAD_NUM];
  Bi bi[THREAD_NUM];
  ui32 pre_p = 0;
  for (ui8 i = 0; i < THREAD_NUM; ++i) {
    ui32 p = (i + 1) * length / THREAD_NUM;
    while (p > 0 && g_buf[p - 1] != '\n') --p;
    bi[i].l = pre_p;
    bi[i].r = p;
    pre_p = p;
    bi[i].tid = i;
    threads[i] = thread(ReadTestData, ref(bi[i]));
  }
  for (auto& thread : threads) {
    thread.join();
  }
#ifdef TEST
  auto end = std::chrono::steady_clock::now();
  double dr_ms = std::chrono::duration<double, std::milli>(end - start).count();
  cout << "[read test data time]: " << dr_ms << " ms\n";
#endif
}

double g_avg_d;
ui32 g_max_edge;
const ui32 g_max_node_num = 2500000;
const ui32 g_max_edge_num = 2500000;
ui16 g_ind[g_max_node_num];
ui16 g_outd[g_max_node_num];
ui32 g_ind_bg[g_max_node_num];
ui32 g_head[g_max_node_num];
struct Edge {
  ui32 to;
  ui32 money;
};
Edge g_edge[g_max_edge_num];

void BuildGraph() {
#ifdef TEST
  auto start = std::chrono::steady_clock::now();
#endif
  vector<ui32> node_d;
  for (auto& it : g_node_data) {
    node_d.insert(node_d.end(), it.begin(), it.end());
  }
  vector<ui32> tmp = node_d;
  sort(tmp.begin(), tmp.end());
  tmp.erase(unique(tmp.begin(), tmp.end()), tmp.end());
  g_node_cnt = tmp.size();
  g_unhash_id = tmp;
  unordered_map<ui32, ui32> hash_id;
  for (ui32 i = 0; i < g_node_cnt; ++i) hash_id[tmp[i]] = i;
  ui32* curlen = new ui32[g_node_cnt + 1]();
  for (auto& it : g_node_data) {
    for (ui32 j = 0; j < it.size(); j += 2) {
      g_outd[hash_id[it[j]]]++;
    }
  }
  g_head[0] = 0;
  for (ui32 i = 1; i <= g_node_cnt; ++i) {
    g_head[i] = g_head[i - 1] + g_outd[i - 1];
  }
  ui32 from, to, money;
  for (ui8 i = 0; i < THREAD_NUM; ++i) {
    for (ui32 j = 0, k = 0; j < g_node_data[i].size(); j += 2, ++k) {
      from = hash_id[g_node_data[i][j]];
      to = hash_id[g_node_data[i][j + 1]];
      money = g_edge_data[i][k];
      g_edge[g_head[from] + curlen[from]].to = to;
      g_edge[g_head[from] + curlen[from]++].money = money;
      g_ind[to]++;
      g_max_edge = g_max_edge < money ? money : g_max_edge;
    }
  }
  g_ind_bg[0] = 0;
  for (ui32 i = 1; i < g_node_cnt; ++i) {
    g_ind_bg[i] = g_ind_bg[i - 1] + g_ind[i - 1];
  }
  g_avg_d =
      (g_ind_bg[g_node_cnt - 1] + g_ind[g_node_cnt - 1]) * 1.0 / g_node_cnt;
#ifdef TEST
  auto end = std::chrono::steady_clock::now();
  double dr_ms = std::chrono::duration<double, std::milli>(end - start).count();
  cout << "[build graph time]: " << dr_ms << " ms\n";
#endif
}

#ifdef TEST
atomic<ui32> cnt(0);
#endif
atomic<ui32> g_th_node(0);
double* g_bc[THREAD_NUM];

class Use16ForDst {
 public:
  struct Heap {
    ui32 num;
    ui16 dist;
    bool operator<(const Heap& t) const { return dist > t.dist; }
    Heap() : num(), dist() {}
    Heap(const ui32& a, const ui16& b) : num(a), dist(b) {}
  };

 public:
  static void DijkWithHeap(const ui32& s, const ui8& tid, ui16* dst,
                           ui16* sigma, double* delta, ui32* lst, bool* visit,
                           Heap* heap) {
    ui32 len_of_lst = 0;
    sigma[s] = 1;
    heap[0] = Heap(s, 0);
    ui32 heap_sz = 1;
    while (heap_sz) {
      pop_heap(heap, heap + heap_sz--);
      if (visit[heap[heap_sz].num]) continue;
      Heap u = heap[heap_sz];
      const ui32& v = u.num;
      visit[v] = true;
      lst[len_of_lst++] = v;
      const ui32 &l = g_head[v], &r = g_head[v + 1];
      Edge* e = &g_edge[l];
      for (ui32 k = l; k < r; ++k, ++e) {
        const ui32& adj_node = e->to;
        ui16 new_dist = u.dist + e->money;
        if (new_dist > dst[adj_node]) continue;
        if (new_dist == dst[adj_node]) {
          sigma[adj_node] += sigma[v];
        } else if (!visit[adj_node]) {
          dst[adj_node] = new_dist;
          heap[heap_sz++] = Heap(adj_node, new_dist);
          push_heap(heap, heap + heap_sz);
          sigma[adj_node] = sigma[v];
        }
      }
    }

    for (ui32 i = len_of_lst - 1; i > 0; --i) {
      const ui32& u = lst[i];
      const ui32 &l = g_head[u], &r = g_head[u + 1];
      Edge* e = &g_edge[l];
      for (ui32 j = l; j < r; ++j, ++e) {
        const ui32& v = e->to;
        if (e->money + dst[u] == dst[v]) {
          delta[u] += (1.0 + delta[v]) * sigma[u] / sigma[v];
        }
      }
      g_bc[tid][u] += delta[u];
    }

    for (ui32 i = 0; i < len_of_lst; ++i) {
      const ui32& v = lst[i];
      dst[v] = INT16_MAX;
      sigma[v] = 0;
      delta[v] = 0.0;
      visit[v] = false;
    }
  }

  static void GetBc(const ui8 tid) {
    bool* visit = new bool[g_node_cnt]();
    ui16* sigma = new ui16[g_node_cnt]();
    ui32* lst = new ui32[g_node_cnt];
    ui16* dst = new ui16[g_node_cnt];
    auto* delta = new double[g_node_cnt]();
    Heap* heap = new Heap[g_max_node_num];
    for (ui32 i = 0; i < g_node_cnt; ++i) {
      dst[i] = INT16_MAX;
    }
    while (true) {
      ui32 node = g_th_node++;
      if (node >= g_node_cnt) break;
      DijkWithHeap(node, tid, dst, sigma, delta, lst, visit, heap);
#ifdef TEST
      ui32 t = ++cnt;
      if (t % 1000 == 0 || t == g_node_cnt) printf("%d/%d\n", t, g_node_cnt);
#endif
    }
    delete[] visit;
    delete[] sigma;
    delete[] lst;
    delete[] dst;
    delete[] delta;
    delete[] heap;
  }

  static void AllocTask() {
#ifdef TEST
    auto start = std::chrono::steady_clock::now();
#endif
    thread threads[THREAD_NUM];
    for (ui8 i = 0; i < THREAD_NUM; ++i) {
      g_bc[i] = new double[g_node_cnt]();
      threads[i] = thread(Use16ForDst::GetBc, i);
    }
    for (auto& thread : threads) {
      thread.join();
    }
#ifdef TEST
    auto end = std::chrono::steady_clock::now();
    double dr_ms =
        std::chrono::duration<double, std::milli>(end - start).count();
    cout << "[calculate time]: " << dr_ms << " ms\n";
#endif
  }
};

class Use16ForDstDense {
 public:
  struct Heap {
    ui32 num;
    ui16 dist;
    bool operator<(const Heap& t) const { return dist > t.dist; }
    Heap() : num(), dist() {}
    Heap(const ui32& a, const ui16& b) : num(a), dist(b) {}
  };

 public:
  static void DijkWithHeap(const ui32& s, const ui8& tid, ui32* pred,
                           ui16* len_of_pred, ui16* dst, ui16* sigma,
                           double* delta, ui32* lst, bool* visit, Heap* heap) {
    ui32 len_of_lst = 0;
    sigma[s] = 1;
    heap[0] = Heap(s, 0);
    ui32 heap_sz = 1;
    while (heap_sz) {
      pop_heap(heap, heap + heap_sz--);
      if (visit[heap[heap_sz].num]) continue;
      Heap u = heap[heap_sz];
      const ui32& v = u.num;
      visit[v] = true;
      lst[len_of_lst++] = v;
      const ui32 &l = g_head[v], &r = g_head[v + 1];
      Edge* e = &g_edge[l];
      for (ui32 k = l; k < r; ++k, ++e) {
        const ui32& adj_node = e->to;
        ui16 new_dist = u.dist + e->money;
        if (new_dist > dst[adj_node]) continue;
        if (new_dist < dst[adj_node]) {
          if (!visit[adj_node]) {
            dst[adj_node] = new_dist;
            heap[heap_sz++] = Heap(adj_node, new_dist);
            push_heap(heap, heap + heap_sz);
            len_of_pred[adj_node] = 0;
            pred[g_ind_bg[adj_node] + len_of_pred[adj_node]++] = v;
            sigma[adj_node] = sigma[v];
          }
        } else {
          pred[g_ind_bg[adj_node] + len_of_pred[adj_node]++] = v;
          sigma[adj_node] += sigma[v];
        }
      }
    }

    for (ui32 i = len_of_lst - 1; i > 0; --i) {
      const ui32& v = lst[i];
      for (ui32 j = 0; j < len_of_pred[v]; ++j) {
        delta[pred[g_ind_bg[v] + j]] +=
            (1.0 + delta[v]) * sigma[pred[g_ind_bg[v] + j]] / sigma[v];
      }
      g_bc[tid][v] += delta[v];
    }

    for (ui32 i = 0; i < len_of_lst; ++i) {
      const ui32& v = lst[i];
      dst[v] = INT16_MAX;
      sigma[v] = 0;
      delta[v] = 0.0;
      visit[v] = false;
    }
  }

  static void GetBc(const ui8 tid) {
    bool* visit = new bool[g_node_cnt]();
    ui16* sigma = new ui16[g_node_cnt]();
    ui32* lst = new ui32[g_node_cnt];
    ui16* dst = new ui16[g_node_cnt];
    ui16* len_of_pred = new ui16[g_node_cnt];
    ui32* pred = new ui32[g_max_edge_num];
    auto* delta = new double[g_node_cnt]();
    Heap* heap = new Heap[g_max_node_num];
    for (ui32 i = 0; i < g_node_cnt; ++i) {
      dst[i] = INT16_MAX;
    }
    while (true) {
      ui32 node = g_th_node++;
      if (node >= g_node_cnt) break;
      DijkWithHeap(node, tid, pred, len_of_pred, dst, sigma, delta, lst, visit,
                   heap);
#ifdef TEST
      ui32 t = ++cnt;
      if (t % 1000 == 0 || t == g_node_cnt) printf("%d/%d\n", t, g_node_cnt);
#endif
    }
    delete[] visit;
    delete[] sigma;
    delete[] lst;
    delete[] dst;
    delete[] pred;
    delete[] len_of_pred;
    delete[] delta;
    delete[] heap;
  }

  static void AllocTask() {
#ifdef TEST
    auto start = std::chrono::steady_clock::now();
#endif
    thread threads[THREAD_NUM];
    for (ui8 i = 0; i < THREAD_NUM; ++i) {
      g_bc[i] = new double[g_node_cnt]();
      threads[i] = thread(Use16ForDstDense::GetBc, i);
    }
    for (auto& thread : threads) {
      thread.join();
    }
#ifdef TEST
    auto end = std::chrono::steady_clock::now();
    double dr_ms =
        std::chrono::duration<double, std::milli>(end - start).count();
    cout << "[calculate time]: " << dr_ms << " ms\n";
#endif
  }
};

class Use32ForDst {
 public:
  struct Heap {
    ui32 num;
    ui32 dist;
    bool operator<(const Heap& t) const { return dist > t.dist; }
    Heap() : num(), dist() {}
    Heap(const ui32& a, const ui32& b) : num(a), dist(b) {}
  };

 public:
  static void DijkWithHeap(const ui32& s, const ui8& tid, ui32* dst,
                           ui16* sigma, double* delta, ui32* lst, bool* visit,
                           Heap* heap) {
    ui32 len_of_lst = 0;
    sigma[s] = 1;
    heap[0] = Heap(s, 0);
    ui32 heap_sz = 1;
    while (heap_sz) {
      pop_heap(heap, heap + heap_sz--);
      if (visit[heap[heap_sz].num]) continue;
      Heap u = heap[heap_sz];
      const ui32& v = u.num;
      visit[v] = true;
      lst[len_of_lst++] = v;
      const ui32 &l = g_head[v], &r = g_head[v + 1];
      Edge* e = &g_edge[l];
      for (ui32 k = l; k < r; ++k, ++e) {
        const ui32& adj_node = e->to;
        ui32 new_dist = u.dist + e->money;
        if (new_dist > dst[adj_node]) continue;
        if (new_dist == dst[adj_node]) {
          sigma[adj_node] += sigma[v];
        } else if (!visit[adj_node]) {
          dst[adj_node] = new_dist;
          heap[heap_sz++] = Heap(adj_node, new_dist);
          push_heap(heap, heap + heap_sz);
          sigma[adj_node] = sigma[v];
        }
      }
    }

    for (ui32 i = len_of_lst - 1; i > 0; --i) {
      const ui32& u = lst[i];
      const ui32 &l = g_head[u], &r = g_head[u + 1];
      Edge* e = &g_edge[l];
      for (ui32 j = l; j < r; ++j, ++e) {
        const ui32& v = e->to;
        if (e->money + dst[u] == dst[v]) {
          delta[u] += (1.0 + delta[v]) * sigma[u] / sigma[v];
        }
      }
      g_bc[tid][u] += delta[u];
    }

    for (ui32 i = 0; i < len_of_lst; ++i) {
      const ui32& v = lst[i];
      dst[v] = INT32_MAX;
      sigma[v] = 0;
      delta[v] = 0.0;
      visit[v] = false;
    }
  }

  static void GetBc(const ui8 tid) {
    bool* visit = new bool[g_node_cnt]();
    ui16* sigma = new ui16[g_node_cnt]();
    ui32* lst = new ui32[g_node_cnt];
    ui32* dst = new ui32[g_node_cnt];
    auto* delta = new double[g_node_cnt]();
    Heap* heap = new Heap[g_max_node_num];
    for (ui32 i = 0; i < g_node_cnt; ++i) {
      dst[i] = INT32_MAX;
    }
    while (true) {
      ui32 node = g_th_node++;
      if (node >= g_node_cnt) break;
      DijkWithHeap(node, tid, dst, sigma, delta, lst, visit, heap);
#ifdef TEST
      ui32 t = ++cnt;
      if (t % 1000 == 0 || t == g_node_cnt) printf("%d/%d\n", t, g_node_cnt);
#endif
    }
    delete[] visit;
    delete[] sigma;
    delete[] lst;
    delete[] dst;
    delete[] delta;
    delete[] heap;
  }

  static void AllocTask() {
#ifdef TEST
    auto start = std::chrono::steady_clock::now();
#endif
    thread threads[THREAD_NUM];
    for (ui8 i = 0; i < THREAD_NUM; ++i) {
      g_bc[i] = new double[g_node_cnt]();
      threads[i] = thread(Use32ForDst::GetBc, i);
    }
    for (auto& thread : threads) {
      thread.join();
    }
#ifdef TEST
    auto end = std::chrono::steady_clock::now();
    double dr_ms =
        std::chrono::duration<double, std::milli>(end - start).count();
    cout << "[calculate time]: " << dr_ms << " ms\n";
#endif
  }
};

class Use32ForDstDense {
 public:
  struct Heap {
    ui32 num;
    ui32 dist;
    bool operator<(const Heap& t) const { return dist > t.dist; }
    Heap() : num(), dist() {}
    Heap(const ui32& a, const ui32& b) : num(a), dist(b) {}
  };

 public:
  static void DijkWithHeap(const ui32& s, const ui8& tid, ui32* pred,
                           ui16* len_of_pred, ui32* dst, ui16* sigma,
                           double* delta, ui32* lst, bool* visit, Heap* heap) {
    ui32 len_of_lst = 0;
    sigma[s] = 1;
    heap[0] = Heap(s, 0);
    ui32 heap_sz = 1;
    while (heap_sz) {
      pop_heap(heap, heap + heap_sz--);
      if (visit[heap[heap_sz].num]) continue;
      Heap u = heap[heap_sz];
      const ui32& v = u.num;
      visit[v] = true;
      lst[len_of_lst++] = v;
      const ui32 &l = g_head[v], &r = g_head[v + 1];
      Edge* e = &g_edge[l];
      for (ui32 k = l; k < r; ++k, ++e) {
        const ui32& adj_node = e->to;
        ui32 new_dist = u.dist + e->money;
        if (new_dist > dst[adj_node]) continue;
        if (new_dist < dst[adj_node]) {
          if (!visit[adj_node]) {
            dst[adj_node] = new_dist;
            heap[heap_sz++] = Heap(adj_node, new_dist);
            push_heap(heap, heap + heap_sz);
            len_of_pred[adj_node] = 0;
            pred[g_ind_bg[adj_node] + len_of_pred[adj_node]++] = v;
            sigma[adj_node] = sigma[v];
          }
        } else {
          pred[g_ind_bg[adj_node] + len_of_pred[adj_node]++] = v;
          sigma[adj_node] += sigma[v];
        }
      }
    }

    for (ui32 i = len_of_lst - 1; i > 0; --i) {
      const ui32& v = lst[i];
      for (ui32 j = 0; j < len_of_pred[v]; ++j) {
        delta[pred[g_ind_bg[v] + j]] +=
            (1.0 + delta[v]) * sigma[pred[g_ind_bg[v] + j]] / sigma[v];
      }
      g_bc[tid][v] += delta[v];
    }

    for (ui32 i = 0; i < len_of_lst; ++i) {
      const ui32& v = lst[i];
      dst[v] = INT32_MAX;
      sigma[v] = 0;
      delta[v] = 0.0;
      visit[v] = false;
    }
  }

  static void GetBc(const ui8 tid) {
    bool* visit = new bool[g_node_cnt]();
    ui16* sigma = new ui16[g_node_cnt]();
    ui32* lst = new ui32[g_node_cnt];
    ui32* dst = new ui32[g_node_cnt];
    ui16* len_of_pred = new ui16[g_node_cnt];
    ui32* pred = new ui32[g_max_edge_num];
    auto* delta = new double[g_node_cnt]();
    Heap* heap = new Heap[g_max_node_num];
    for (ui32 i = 0; i < g_node_cnt; ++i) {
      dst[i] = INT32_MAX;
    }
    while (true) {
      ui32 node = g_th_node++;
      if (node >= g_node_cnt) break;
      DijkWithHeap(node, tid, pred, len_of_pred, dst, sigma, delta, lst, visit,
                   heap);
#ifdef TEST
      ui32 t = ++cnt;
      if (t % 1000 == 0 || t == g_node_cnt) printf("%d/%d\n", t, g_node_cnt);
#endif
    }
    delete[] visit;
    delete[] sigma;
    delete[] lst;
    delete[] dst;
    delete[] pred;
    delete[] len_of_pred;
    delete[] delta;
    delete[] heap;
  }

  static void AllocTask() {
#ifdef TEST
    auto start = std::chrono::steady_clock::now();
#endif
    thread threads[THREAD_NUM];
    for (ui8 i = 0; i < THREAD_NUM; ++i) {
      g_bc[i] = new double[g_node_cnt]();
      threads[i] = thread(Use32ForDstDense::GetBc, i);
    }
    for (auto& thread : threads) {
      thread.join();
    }
#ifdef TEST
    auto end = std::chrono::steady_clock::now();
    double dr_ms =
        std::chrono::duration<double, std::milli>(end - start).count();
    cout << "[calculate time]: " << dr_ms << " ms\n";
#endif
  }
};

class Use64ForDst {
 public:
  struct Heap {
    ui32 num;
    ui64 dist;
    bool operator<(const Heap& t) const { return dist > t.dist; }
    Heap() : num(), dist() {}
    Heap(const ui32& a, const ui64& b) : num(a), dist(b) {}
  };

 public:
  static void DijkWithHeap(const ui32& s, const ui8& tid, ui64* dst,
                           ui16* sigma, double* delta, ui32* lst, bool* visit,
                           Heap* heap) {
    ui32 len_of_lst = 0;
    sigma[s] = 1;
    heap[0] = Heap(s, 0);
    ui32 heap_sz = 1;
    while (heap_sz) {
      pop_heap(heap, heap + heap_sz--);
      if (visit[heap[heap_sz].num]) continue;
      Heap u = heap[heap_sz];
      const ui32& v = u.num;
      visit[v] = true;
      lst[len_of_lst++] = v;
      const ui32 &l = g_head[v], &r = g_head[v + 1];
      Edge* e = &g_edge[l];
      for (ui32 k = l; k < r; ++k, ++e) {
        const ui32& adj_node = e->to;
        ui64 new_dist = u.dist + e->money;
        if (new_dist > dst[adj_node]) continue;
        if (new_dist == dst[adj_node]) {
          sigma[adj_node] += sigma[v];
        } else if (!visit[adj_node]) {
          dst[adj_node] = new_dist;
          heap[heap_sz++] = Heap(adj_node, new_dist);
          push_heap(heap, heap + heap_sz);
          sigma[adj_node] = sigma[v];
        }
      }
    }

    for (ui32 i = len_of_lst - 1; i > 0; --i) {
      const ui32& u = lst[i];
      const ui32 &l = g_head[u], &r = g_head[u + 1];
      Edge* e = &g_edge[l];
      for (ui32 j = l; j < r; ++j, ++e) {
        const ui32& v = e->to;
        if (e->money + dst[u] == dst[v]) {
          delta[u] += (1.0 + delta[v]) * sigma[u] / sigma[v];
        }
      }
      g_bc[tid][u] += delta[u];
    }

    for (ui32 i = 0; i < len_of_lst; ++i) {
      const ui32& v = lst[i];
      dst[v] = INT64_MAX;
      sigma[v] = 0;
      delta[v] = 0.0;
      visit[v] = false;
    }
  }

  static void GetBc(const ui8 tid) {
    bool* visit = new bool[g_node_cnt]();
    ui16* sigma = new ui16[g_node_cnt]();
    ui32* lst = new ui32[g_node_cnt];
    ui64* dst = new ui64[g_node_cnt];
    auto* delta = new double[g_node_cnt]();
    Heap* heap = new Heap[g_max_node_num];
    for (ui32 i = 0; i < g_node_cnt; ++i) {
      dst[i] = INT64_MAX;
    }
    while (true) {
      ui32 node = g_th_node++;
      if (node >= g_node_cnt) break;
      DijkWithHeap(node, tid, dst, sigma, delta, lst, visit, heap);
#ifdef TEST
      ui32 t = ++cnt;
      if (t % 1000 == 0 || t == g_node_cnt) printf("%d/%d\n", t, g_node_cnt);
#endif
    }
    delete[] visit;
    delete[] sigma;
    delete[] lst;
    delete[] dst;
    delete[] delta;
    delete[] heap;
  }

  static void AllocTask() {
#ifdef TEST
    auto start = std::chrono::steady_clock::now();
#endif
    thread threads[THREAD_NUM];
    for (ui8 i = 0; i < THREAD_NUM; ++i) {
      g_bc[i] = new double[g_node_cnt]();
      threads[i] = thread(Use64ForDst::GetBc, i);
    }
    for (auto& thread : threads) {
      thread.join();
    }
#ifdef TEST
    auto end = std::chrono::steady_clock::now();
    double dr_ms =
        std::chrono::duration<double, std::milli>(end - start).count();
    cout << "[calculate time]: " << dr_ms << " ms\n";
#endif
  }
};

class Use64ForDstDense {
 public:
  struct Heap {
    ui32 num;
    ui64 dist;
    bool operator<(const Heap& t) const { return dist > t.dist; }
    Heap() : num(), dist() {}
    Heap(const ui32& a, const ui64& b) : num(a), dist(b) {}
  };

 public:
  static void DijkWithHeap(const ui32& s, const ui8& tid, ui32* pred,
                           ui16* len_of_pred, ui64* dst, ui16* sigma,
                           double* delta, ui32* lst, bool* visit, Heap* heap) {
    ui32 len_of_lst = 0;
    sigma[s] = 1;
    heap[0] = Heap(s, 0);
    ui32 heap_sz = 1;
    while (heap_sz) {
      pop_heap(heap, heap + heap_sz--);
      if (visit[heap[heap_sz].num]) continue;
      Heap u = heap[heap_sz];
      const ui32& v = u.num;
      visit[v] = true;
      lst[len_of_lst++] = v;
      const ui32 &l = g_head[v], &r = g_head[v + 1];
      Edge* e = &g_edge[l];
      for (ui32 k = l; k < r; ++k, ++e) {
        const ui32& adj_node = e->to;
        ui64 new_dist = u.dist + e->money;
        if (new_dist > dst[adj_node]) continue;
        if (new_dist < dst[adj_node]) {
          if (!visit[adj_node]) {
            dst[adj_node] = new_dist;
            heap[heap_sz++] = Heap(adj_node, new_dist);
            push_heap(heap, heap + heap_sz);
            len_of_pred[adj_node] = 0;
            pred[g_ind_bg[adj_node] + len_of_pred[adj_node]++] = v;
            sigma[adj_node] = sigma[v];
          }
        } else {
          pred[g_ind_bg[adj_node] + len_of_pred[adj_node]++] = v;
          sigma[adj_node] += sigma[v];
        }
      }
    }

    for (ui32 i = len_of_lst - 1; i > 0; --i) {
      const ui32& v = lst[i];
      for (ui32 j = 0; j < len_of_pred[v]; ++j) {
        delta[pred[g_ind_bg[v] + j]] +=
            (1.0 + delta[v]) * sigma[pred[g_ind_bg[v] + j]] / sigma[v];
      }
      g_bc[tid][v] += delta[v];
    }

    for (ui32 i = 0; i < len_of_lst; ++i) {
      const ui32& v = lst[i];
      dst[v] = INT64_MAX;
      sigma[v] = 0;
      delta[v] = 0.0;
      visit[v] = false;
    }
  }

  static void GetBc(const ui8 tid) {
    bool* visit = new bool[g_node_cnt]();
    ui16* sigma = new ui16[g_node_cnt]();
    ui32* lst = new ui32[g_node_cnt];
    ui64* dst = new ui64[g_node_cnt];
    ui16* len_of_pred = new ui16[g_node_cnt];
    ui32* pred = new ui32[g_max_edge_num];
    auto* delta = new double[g_node_cnt]();
    Heap* heap = new Heap[g_max_node_num];
    for (ui32 i = 0; i < g_node_cnt; ++i) {
      dst[i] = INT64_MAX;
    }
    while (true) {
      ui32 node = g_th_node++;
      if (node >= g_node_cnt) break;
      DijkWithHeap(node, tid, pred, len_of_pred, dst, sigma, delta, lst, visit,
                   heap);
#ifdef TEST
      ui32 t = ++cnt;
      if (t % 1000 == 0 || t == g_node_cnt) printf("%d/%d\n", t, g_node_cnt);
#endif
    }
    delete[] visit;
    delete[] sigma;
    delete[] lst;
    delete[] dst;
    delete[] pred;
    delete[] len_of_pred;
    delete[] delta;
    delete[] heap;
  }

  static void AllocTask() {
#ifdef TEST
    auto start = std::chrono::steady_clock::now();
#endif
    thread threads[THREAD_NUM];
    for (ui8 i = 0; i < THREAD_NUM; ++i) {
      g_bc[i] = new double[g_node_cnt]();
      threads[i] = thread(Use64ForDstDense::GetBc, i);
    }
    for (auto& thread : threads) {
      thread.join();
    }
#ifdef TEST
    auto end = std::chrono::steady_clock::now();
    double dr_ms =
        std::chrono::duration<double, std::milli>(end - start).count();
    cout << "[calculate time]: " << dr_ms << " ms\n";
#endif
  }
};

Inline bool MyCmp(const pair<double, ui32>& a, const pair<double, ui32>& b) {
  if (abs(a.first - b.first) > 0.0001)
    return a.first > b.first;
  else
    return a.second < b.second;
}

int main() {
  const string test_file = "/data/test_data.txt";
  const string result_file = "/projects/student/result.txt";
#ifdef TEST
  auto start = std::chrono::steady_clock::now();
#endif
  Read(test_file);
  BuildGraph();
  if (g_max_edge < 256) {
    if (g_avg_d < 2) {
      Use16ForDst::AllocTask();
    } else {
      Use16ForDstDense::AllocTask();
    }
  } else if (g_max_edge < 1000001) {
    if (g_avg_d < 2) {
      Use32ForDst::AllocTask();
    } else {
      Use32ForDstDense::AllocTask();
    }
  } else {
    if (g_avg_d < 2) {
      Use64ForDst::AllocTask();
    } else {
      Use64ForDstDense::AllocTask();
    }
  }
  auto* bc_final = new double[g_node_cnt]();
  for (auto& ii : g_bc) {
    for (ui32 j = 0; j < g_node_cnt; ++j) {
      bc_final[j] += ii[j];
    }
  }
  auto res = new pair<double, ui32>[g_node_cnt];
  for (ui32 i = 0; i < g_node_cnt; ++i) {
    res[i] = make_pair(bc_final[i], i);
  }
  sort(res, res + g_node_cnt, MyCmp);
  ui8 index = g_node_cnt < 100 ? g_node_cnt : 100;
  ofstream write(result_file.c_str(), ios::out);
  for (ui8 i = 0; i < index; ++i) {
    write << g_unhash_id[res[i].second] << ',' << fixed << setprecision(3)
          << res[i].first << '\n';
  }
  write.close();
  delete[] bc_final;
  delete[] res;
#ifdef TEST
  auto end = std::chrono::steady_clock::now();
  double dr_ms = std::chrono::duration<double, std::milli>(end - start).count();
  cout << "[total cost time]: " << dr_ms << " ms\n";
#endif
  return 0;
}
