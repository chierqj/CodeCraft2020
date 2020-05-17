//  Copyright (C) 2020 Mingliang Zeng (mlzeng.com)

//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.

//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.

//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <bitset>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <map>
#include <random>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

// options for tuning performance
#define USE_MMAP_INPUT
#define ALLOW_UINT16
#define USE_LSH
#define ALLOW_TOPO_PRUNE

// implies
#ifdef ALLOW_TOPO_PRUNE
#define USE_LSH
#endif

// options for debug only
#define INFO
#undef STATISTICS

// additional functions
#define WEIGHT
#undef KASHI

// parameters
#ifdef LOCAL_TEST
const bool submit = false;
#else
const bool submit = true;
#endif
const int num_threads = 4;  // affect total_size
const int num_parts = 1024;
const int maxe = 2100000;  // max edges
const int maxn = 2100000;  // max nodes
const int maxb = 1000000;  // reverse answers buffer size
const int minl = 3;
const int maxl = 7;
const int avgd = 20;  // average deg
const int int_to_str_align = 16;
const int len_pos = 11;
const int basket_size = 4096;  // input basket
const int factor_A = 3;
const int factor_B = 5;
const char *input_file_name =
    (submit ? "/data/test_data.txt" : "test_data.txt");
const char *output_file_name =
    (submit ? "/projects/student/result.txt" : "output.txt");
const size_t block_size = ((size_t)1) << (30);
const size_t total_size = block_size << (2 + 3);  // tid:2 len:3
const int pre_sleep_time = 0;
const int post_sleep_time = 0;
const int fixed_len[] = {0, 12, 24, 36, 48, 56, 72, 80};
// 0:disable  1:simple  2:full
const std::vector<int> topo_prune_modes_list = {2};
// 00:xx  01:xA  10:Bx  11:BA
const std::vector<int> adj_list_2_modes_list = {3};
std::random_device my_random_device;
std::default_random_engine rng(my_random_device());
const int topo_prune_mode =
    topo_prune_modes_list[std::uniform_int_distribution<int>(
        0, topo_prune_modes_list.size() - 1)(rng)];
const int adj_list_2_mode =
    adj_list_2_modes_list[std::uniform_int_distribution<int>(
        0, adj_list_2_modes_list.size() - 1)(rng)];

typedef uint8_t ts_type;  // uint8_t or uint16_t

#define clock_type steady_clock
#define likely(x) (__builtin_expect(!!(x), 1))
#define unlikely(x) (__builtin_expect(!!(x), 0))

enum SP {
  SP_INPUT,
  SP_LSH,
  SP_PSRS,
  SP_UNIQ,
  SP_PRUNE,
  SP_DEG,
  SP_REM,
  SP_RLB,
  SP_DS,
  SP_DS_0,
  SP_DS_1,
  SP_DS_2,
  SP_DS_3,
  SP_DS_4,
  SP_DS_5,
  SP_DS_6,
  SP_DS_7,
  SP_RUN,
  SP_STEP1,
  SP_STEP2,
  SP_UPD,
  SP_OUTPUT,
  SP_TOTAL,
  SP_MAX = SP_TOTAL
};

template <typename T>
struct Edge {
  T to;
#ifdef WEIGHT
  int w1;
#endif
  inline const bool operator<(const Edge &rhs) const { return to < rhs.to; }
  inline const bool operator==(const Edge &rhs) const { return to == rhs.to; }
};

template <typename T>
struct Edge2 : Edge<T> {
  T im;
#ifdef WEIGHT
  int w2;
#endif
};

struct Node {
  int L, L2, R, R2, C, C2;
};

struct InputItem {
  uint32_t ofs, val;
  inline const bool operator<(const InputItem &rhs) const {
    return val < rhs.val;
  }
};

template <typename T>
struct Mark {
  T s1;
#ifdef WEIGHT
  int w1;
#endif
  T s2;
  T a2;
#ifdef WEIGHT
  int w2s, w2t;
#endif
  bool c2;
  int l2, r2;
  T s3;
  T a3;
  T b3;
#ifdef WEIGHT
  int w3s, w3t;
#endif
  bool c3;
  int l3, r3;
  inline void update1(T t, int w) {
    s1 = t;
#ifdef WEIGHT
    w1 = w;
#endif
  }
  inline void update2(T t, T a, int ws, int wt) {
    if (likely(s2 != t)) {
      s2 = t;
      a2 = a;
      c2 = false;
#ifdef WEIGHT
      w2s = ws;
      w2t = wt;
#endif
    } else {
      r2 = 0;
      c2 = true;
    }
  }
  inline void update3(T t, T a, T b, int ws, int wt) {
    if (likely(s3 != t)) {
      s3 = t;
      a3 = a;
      b3 = b;
      c3 = false;
#ifdef WEIGHT
      w3s = ws;
      w3t = wt;
#endif
    } else {
      r3 = 0;
      c3 = true;
    }
  }
  inline bool check1(T t) { return s1 == t; }
  inline bool check2(T t) { return s2 == t; }
  inline bool check3(T t) { return s3 == t; }
};

template <typename T, uint32_t m>
inline uint32_t my_hash(const T &a) {
  return (uint32_t)std::hash<T>{}(a) % m;
}

// global variables
static char *__restrict__ region;
static bool finished[maxn];
static Edge<uint32_t> _LA[maxe];
static Edge<uint32_t> _LB[maxe];
static Edge2<uint32_t> _LA2[maxe * avgd];
static Edge2<uint32_t> _LB2[maxe * avgd];
static Node NA[maxn];
static Node NB[maxn];
static char int_to_str[maxn][int_to_str_align];
static int nodes_count;
static int edges_count;
static std::vector<int> edges_count_each;
static int edges_list[maxe][4];
static InputItem ii_list[maxe * 2];
static int id_rev[maxe * 2];
static int in_fd, out_fd;
static char *answers;
static char *answers_ptr_buf[num_threads][maxl + 1];
static std::pair<char *, char *> answers_ptr[num_parts][maxl + 1];
static size_t input_file_size;
static size_t output_text_size;
static Mark<uint32_t> _marks[num_threads][maxn];
static ts_type _ts[num_threads][maxn];
static uint32_t _VA[maxn];
static uint32_t _VB[maxn];
static uint32_t _sorted[maxe * 2];
static uint32_t _mark_answers[num_threads][maxb];
static size_t answers_text_start_point[num_parts][maxl + 1];
static std::atomic<int> answers_count;

const char cm = ',';
const char rt = '\r';
const char nl = '\n';

// probability of having valid weights (paths of length n)
const double path_prob[] = {1,        1,        0.733333, 0.576395,
                            0.456018, 0.361079, 0.285957, 0.226473};

// probability of having valid weights (cycles of length n)
const double cycle_prob[] = {1,        1,        0.666667, 0.506667,
                             0.395281, 0.312003, 0.246877, 0.191461};

#ifdef STATISTICS
std::atomic<uint64_t> search_count;
std::atomic<uint64_t> hit_count;
std::atomic<uint64_t> fb_count;
#endif

namespace Color {
inline void reset() { fprintf(stderr, "\033[0m"); }
inline void red() { fprintf(stderr, "\033[1;31m"); }
inline void green() { fprintf(stderr, "\033[1;32m"); }
inline void yellow() { fprintf(stderr, "\033[1;33m"); }
inline void blue() { fprintf(stderr, "\033[1;34m"); }
inline void magenta() { fprintf(stderr, "\033[1;35m"); }
inline void cyan() { fprintf(stderr, "\033[1;36m"); }
inline void orange() { fprintf(stderr, "\033[38;5;214m"); }
inline void newline() { fprintf(stderr, "\n"); }
}  // namespace Color

namespace Timer {
static std::array<std::chrono::_V2::clock_type::time_point, maxn>
    start_point[SP_MAX + 1], end_point[SP_MAX + 1];
static std::atomic<int> start_counter[SP_MAX + 1], end_counter[SP_MAX + 1];

inline void start(SP x) {
#ifdef INFO
  start_point[x][std::atomic_fetch_add_explicit(&start_counter[x], 1,
                                                std::memory_order_relaxed)] =
      std::chrono::clock_type::now();
#endif
}

inline void end(SP x) {
#ifdef INFO
  end_point[x][std::atomic_fetch_add_explicit(&end_counter[x], 1,
                                              std::memory_order_relaxed)] =
      std::chrono::clock_type::now();
#endif
}

inline long long gettime(SP x) {
  long long ns = 0;
  int sz = start_counter[x];
  for (int i = 0; i < sz; i++) {
    auto duration = end_point[x][i] - start_point[x][i];
    ns +=
        std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count();
  }
  return ns;
}

inline void print(SP x, const char *s) {
#ifdef INFO
  long long ns = gettime(x);
  long long ms = round(ns / 1e6);
  long long us = round(ns / 1e3);
  fprintf(stderr, "%6s: %-5lld ms  ( %-8lld us )\n", s, ms, us);
#endif
}
}  // namespace Timer

namespace PSRS {
// Sort an array using in parallel using multiple threads
// Copyright (C) 2016 Amir Baserinia (baserinia.github.io)

template <typename T>
static inline void _sort(size_t sz, T *list) {
  std::sort(list, list + sz);
}

template <typename T>
static void _count(size_t sz, const T *list, const T *range, size_t *bucket) {
  const unsigned np = num_threads;
  for (unsigned n = 0; n < np; n++) bucket[n] = 0;
  for (size_t i = 0; i < sz; i++) {
    unsigned n;
    for (n = 0; n < np - 1; n++) {
      if (list[i] < range[n]) {
        bucket[n]++;
        break;
      }
    }
    if (n == np - 1) bucket[np - 1]++;
  }
}

template <typename T>
static void _reorder(size_t sz, const T *list, const T *range, size_t *map,
                     T *out) {
  const unsigned np = num_threads;
  for (size_t i = 0; i < sz; i++) {
    unsigned n;
    for (n = 0; n < np - 1; n++) {
      if (list[i] < range[n]) {
        out[map[n]++] = list[i];
        break;
      }
    }
    if (n == np - 1) out[map[np - 1]++] = list[i];
  }
}

template <typename T>
void parasort(size_t sz, T *list) {
  const unsigned sf = 256;
  const unsigned np = num_threads;
  const unsigned nSample = np * sf;
  const unsigned nMap = np * np;
  static T range[np];
  static T sample[nSample];
  T *sorted = (T *)_sorted;
  static size_t tmp[nMap];
  static size_t map[nMap];
  static size_t bucket[nMap];
  std::array<std::thread, np> threads;
  for (unsigned i = 0; i < nSample; i++)
    sample[i] = list[std::uniform_int_distribution<int>(0, sz - 1)(rng)];
  _sort(nSample, &sample[0]);
  for (unsigned i = 0; i < np - 1; i++) range[i] = sample[(i + 1) * sf];
  unsigned part = sz / np;
  for (unsigned i = 0; i < np; i++) {
    size_t start = i * part;
    size_t length = (i + 1 == np) ? sz - i * part : part;
    threads[i] = std::thread(_count<T>, length, &list[start], &range[0],
                             &bucket[i * np]);
  }
  for (auto &&thread : threads) thread.join();
  for (unsigned i = 0; i < nMap; i++)
    tmp[i] = i ? tmp[i - 1] + bucket[((i - 1) % np) * np + (i - 1) / np] : 0;
  for (unsigned i = 0; i < nMap; i++) map[i] = tmp[(i % np) * np + i / np];
  for (unsigned i = 0; i < nMap; i++) tmp[i] = map[i];
  for (unsigned i = 0; i < np; i++) {
    size_t start = i * part;
    size_t length = (i + 1 == np) ? sz - i * part : part;
    threads[i] = std::thread(_reorder<T>, length, &list[start], &range[0],
                             &tmp[i * np], &sorted[0]);
  }
  for (auto &&thread : threads) thread.join();
  for (unsigned i = 0; i < np; i++) {
    size_t start = map[i];
    size_t length = (i + 1 == np) ? sz - map[i] : map[i + 1] - map[i];
    threads[i] = std::thread(_sort<T>, length, &sorted[start]);
  }
  for (auto &&thread : threads) thread.join();
  for (unsigned i = 0; i < np; i++) {
    size_t start = (sz * (i + 0)) / np;
    size_t end = (sz * (i + 1)) / np;
    threads[i] = std::thread(memcpy, list + start, sorted + start,
                             (end - start) * sizeof(T));
  }
  for (auto &&thread : threads) thread.join();
}
};  // namespace PSRS

inline int get_job() {
  static std::atomic<int> current_job;
  return std::atomic_fetch_add_explicit(&current_job, 1,
                                        std::memory_order_relaxed);
}

template <int nt>
inline void parallel_run(void (*func)(int)) {
  static std::thread threads[nt];
  for (int i = 0; i < nt; i++) {
    threads[i] = std::thread(func, i);
  }
  for (int i = 0; i < nt; i++) {
    threads[i].join();
  }
}

inline bool is_digit(char &c) { return '0' <= c /* && c <= '9'*/; }

inline bool check_weight(int64_t x, int64_t y) {  // prob = 11/15 = 0.733333
  return (factor_B * y >= x) && (factor_A * x >= y);
}

inline char *answers_initial_pointer(int tid, int len) {
  return (answers) + block_size * (tid * 8 + len);
}

void daemon() {
  if (fork()) {
    exit(0);
  }
  setsid();
}

template <int x>
inline bool delay() {
  if (x > 0) {
    std::this_thread::sleep_for(std::chrono::milliseconds(x));
  }
  if (x < 0) {
    return (fork()) ? (false) : (daemon(), true);
  }
  return true;
}

template <typename T>
inline void init_answers() {
  answers =
      (char *)(mmap64(NULL, total_size, PROT_READ | PROT_WRITE,
                      MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0));
  madvise(answers, total_size, MADV_SEQUENTIAL);
}

template <typename T>
inline void update_graph() {
  static std::atomic<int> mtx_graph;
  static int P = 0;
  T *VA = (T *)_VA;
  T *VB = (T *)_VB;
  const Edge<T> *__restrict__ LA = (Edge<T> *)_LA;
  const Edge<T> *__restrict__ LB = (Edge<T> *)_LB;
  Timer::start(SP_UPD);
  if ((mtx_graph++) == 0) {  // try lock
    while (P < nodes_count && finished[P]) {
      int la = 0;
      int lb = 0;
      for (auto e = NA[P].L; e < NA[P].R; e++) {
        VB[lb++] = LA[e].to;
      }
      for (auto e = NB[P].L; e < NB[P].R; e++) {
        VA[la++] = LB[e].to;
      }
      if (adj_list_2_mode & 1) {
        const Edge2<T> *__restrict__ LA2 = (Edge2<T> *)_LA2;
        for (int i = 0; i < la; i++) {
          NA[VA[i]].L++;
          while (LA2[NA[VA[i]].L2].im == P) {
            NA[VA[i]].L2++;
          }
        }
      }
      if (adj_list_2_mode & 2) {
        const Edge2<T> *__restrict__ LB2 = (Edge2<T> *)_LB2;
        for (int i = 0; i < lb; i++) {
          NB[VB[i]].L++;
          while (LB2[NB[VB[i]].L2].im == P) {
            NB[VB[i]].L2++;
          }
        }
      }
      P++;
    }
  }
  mtx_graph--;
  Timer::end(SP_UPD);
}

template <typename T, bool mode>
inline void process_marks(ts_type ts[], Mark<T> marks[], T start) {
  const Edge<T> *__restrict__ LB = (Edge<T> *)_LB;
  if (mode) {
    const Edge2<T> *__restrict__ LB2 = (Edge2<T> *)_LB2;
    for (auto e = NB[start].L; e < NB[start].R; e++) {  // level 1
      auto &u = LB[e].to;
      if (u <= start) {
        continue;
      }
#ifdef WEIGHT
      auto &w1 = LB[e].w1;
#else
      int w1 = 0;
#endif
      ts[u] = start;
      marks[u].update1(start, w1);
    }
    for (auto e = NB[start].L2; e < NB[start].R2; e++) {  // level 2
      auto &u = LB2[e].im;
      auto &v = LB2[e].to;
      if (u <= start || v <= start) {
        continue;
      }
#ifdef WEIGHT
      auto &w1 = LB2[e].w1;
      auto &w2 = LB2[e].w2;
#else
      int w1 = 0;
      int w2 = 0;
#endif
      ts[v] = start;
      marks[v].update2(start, u, w1, w2);
    }
    for (auto e = NB[start].L; e < NB[start].R; e++) {  // level 3
      auto &u = LB[e].to;
      if (u <= start) {
        continue;
      }
#ifdef WEIGHT
      auto &w0 = LB[e].w1;
#else
      int w0 = 0;
#endif
#ifdef WEIGHT
      int64_t max_w1 = factor_B * (int64_t)w0;
      int64_t min_w1 = (factor_A - 1 + (int64_t)w0) / factor_A;
#endif
      for (auto f = NB[u].L2; f < NB[u].R2; f++) {
        auto &v = LB2[f].im;
        auto &w = LB2[f].to;
        if (v <= start || w <= start) {
          continue;
        }
#ifdef WEIGHT
        auto &w1 = LB2[f].w1;
        auto &w2 = LB2[f].w2;
        if (w1 < min_w1 || w1 > max_w1) {
          continue;
        }
#else
        int w1 = 0;
        int w2 = 0;
#endif
        ts[w] = start;
        marks[w].update3(start, v, u, w0, w2);
      }
    }
  } else {
    for (auto e = NB[start].L; e < NB[start].R; e++) {  // level 1
      auto &u = LB[e].to;
      if (u <= start) {
        continue;
      }
#ifdef WEIGHT
      auto &w1 = LB[e].w1;
#else
      int w1 = 0;
#endif
      ts[u] = start;
      marks[u].update1(start, w1);
      for (auto f = NB[u].L; f < NB[u].R; f++) {  // level 2
        auto &v = LB[f].to;
        if (v <= start) {
          continue;
        }
#ifdef WEIGHT
        auto &w2 = LB[f].w1;
        if (!check_weight(w2, w1)) {
          continue;
        }
#else
        int w2 = 0;
#endif
        ts[v] = start;
        marks[v].update2(start, u, w1, w2);
        for (auto g = NB[v].L; g < NB[v].R; g++) {  // level 3
          auto &w = LB[g].to;
          if (w <= start || w == u) {
            continue;
          }
#ifdef WEIGHT
          auto &w3 = LB[g].w1;
          if (!check_weight(w3, w2)) {
            continue;
          }
#else
          int w3 = 0;
#endif
          ts[w] = start;
          marks[w].update3(start, v, u, w1, w3);
        }
      }
    }
  }
}

template <typename T, int d1>
inline void print_answers(Mark<T> marks[], T st[], T start, T u, int pid,
                          char *ans[], int ws, int wt, int &ans_p, int tid,
                          int &ans_cnt, char *line, char *line_ptr) {
  auto mark_answers = _mark_answers[tid];
  auto str_len = line_ptr - line;
  const Edge<T> *__restrict__ LA = (Edge<T> *)_LA;
  const size_t bs = maxl + 1;
#ifdef STATISTICS
  search_count++;
#endif
  auto &m = marks[u];
  {
    if (d1 == 2 && likely(m.check1(start))) {
      const int d2 = 1;
#ifdef WEIGHT
      const bool cw = check_weight(wt, m.w1) && check_weight(m.w1, ws);
#else
      const bool cw = true;
#endif
#ifdef STATISTICS
      hit_count++;
#endif
      if (cw) {  // length = 3
        ans_cnt++;
        memcpy(ans[d1 + d2], line, fixed_len[d1 + 1]);
        ans[d1 + d2] += str_len;
        *(ans[d1 + d2] - 1) = '\n';
      }
    }
    if (likely(m.check2(start))) {  // length = 4 or 6
#ifdef STATISTICS
      hit_count++;
#endif
      const int d2 = 2;
      if (likely(!m.c2)) {
        auto &v = m.a2;
        // check duplication
        bool flag = false;
        for (int i = 1; i <= d1; i++) {
          if (st[i] == v) {
            flag = true;
          }
        }
        if (!flag) {
#ifdef WEIGHT
          const bool cw = check_weight(wt, m.w2t) && check_weight(m.w2s, ws);
#else
          const bool cw = true;
#endif
          if (cw) {
            ans_cnt++;
            memcpy(ans[d1 + d2], line, fixed_len[d1 + 1]);
            ans[d1 + d2] += str_len;
            memcpy(ans[d1 + d2], int_to_str[v], 16);
            ans[d1 + d2] += int_to_str[v][len_pos];
            *(ans[d1 + d2] - 1) = '\n';
          }
        }
      } else {
#ifdef STATISTICS
        fb_count++;
#endif
        if (m.r2 == 0) {
          int o;
          m.l2 = o = ans_p;
          for (auto e = NA[u].L; e < NA[u].R; e++) {
            auto &v = LA[e].to;
            auto &y = marks[v];
            if (y.check1(start)) {
#ifdef WEIGHT
              auto &w1 = LA[e].w1;
              if (!check_weight(w1, y.w1)) {
                continue;
              }
#endif
#ifdef WEIGHT
              mark_answers[o + 0] = w1;
              mark_answers[o + 1] = y.w1;
#endif
              mark_answers[o + 2] = v;
              memcpy(&mark_answers[o + 3], int_to_str[v], 16);
              o += 6;
            }
          }
          m.r2 = ans_p = o;
        }
        {
#ifdef WEIGHT
          int64_t max_w1 = factor_A * (int64_t)wt;
          int64_t min_w1 = (factor_B - 1 + (int64_t)wt) / factor_B;
          int64_t max_w2 = factor_B * (int64_t)ws;
          int64_t min_w2 = (factor_A - 1 + (int64_t)ws) / factor_A;
#endif
          for (int o = m.l2; o < m.r2; o += 6) {
            int w1 = mark_answers[o + 0];
            int w2 = mark_answers[o + 1];
            T v = mark_answers[o + 2];
#ifdef WEIGHT
            if (w1 < min_w1 || w1 > max_w1 || w2 < min_w2 || w2 > max_w2) {
              continue;
            }
#endif
            // check duplication
            bool flag = false;
            for (int i = 1; i <= d1; i++) {
              if (st[i] == v) {
                flag = true;
              }
            }
            if (flag) {
              continue;
            }
            ans_cnt++;
            memcpy(ans[d1 + d2], line, fixed_len[d1 + 1]);
            ans[d1 + d2] += str_len;
            memcpy(ans[d1 + d2], &mark_answers[o + 3], 16);
            ans[d1 + d2] += ((char *)&mark_answers[o + 3])[len_pos];
            *(ans[d1 + d2] - 1) = '\n';
          }
        }
      }
    }
    if (likely(m.check3(start))) {  // length = 5 or 7
#ifdef STATISTICS
      hit_count++;
#endif
      const int d2 = 3;
      if (likely(!m.c3)) {
        auto &v = m.a3;
        auto &w = m.b3;
        // check duplication
        bool flag = false;
        for (int i = 1; i <= d1; i++) {
          if (st[i] == v || st[i] == w) {
            flag = true;
          }
        }
        if (flag) {
          return;
        }
#ifdef WEIGHT
        if (!check_weight(wt, m.w3t) || !check_weight(m.w3s, ws)) {
          return;
        }
#endif
        ans_cnt++;
        memcpy(ans[d1 + d2], line, fixed_len[d1 + 1]);
        ans[d1 + d2] += str_len;
        memcpy(ans[d1 + d2], int_to_str[v], 16);
        ans[d1 + d2] += int_to_str[v][len_pos];
        memcpy(ans[d1 + d2], int_to_str[w], 16);
        ans[d1 + d2] += int_to_str[w][len_pos];
        *(ans[d1 + d2] - 1) = '\n';
      } else {
#ifdef STATISTICS
        fb_count++;
#endif
        if (m.r3 == 0) {
          int o;
          m.l3 = o = ans_p;
          for (auto e = NA[u].L; e < NA[u].R; e++) {
            auto &v = LA[e].to;
            auto &y = marks[v];
            if (y.check2(start)) {
              for (auto f = NA[v].L; f < NA[v].R; f++) {
                auto &w = LA[f].to;
                auto &z = marks[w];
                if (z.check1(start)) {
#ifdef WEIGHT
                  auto &w1 = LA[e].w1;
                  auto &w2 = LA[f].w1;
                  if (!check_weight(w1, w2) || !check_weight(w2, z.w1)) {
                    continue;
                  }
#endif
                  {
#ifdef WEIGHT
                    mark_answers[o + 0] = w1;
                    mark_answers[o + 1] = z.w1;
#endif
                    mark_answers[o + 2] = v;
                    mark_answers[o + 3] = w;
                    memcpy(&mark_answers[o + 4], int_to_str[v],
                           int_to_str_align);
                    memcpy(
                        ((char *)&mark_answers[o + 4]) + int_to_str[v][len_pos],
                        int_to_str[w], int_to_str_align);
                    ((char *)&mark_answers[o + 4])[len_pos * 2] =
                        int_to_str[v][len_pos] + int_to_str[w][len_pos];
                    o += 10;
                  }
                }
              }
            }
          }
          m.r3 = ans_p = o;
        }
        {
#ifdef WEIGHT
          int64_t max_w1 = factor_A * (int64_t)wt;
          int64_t min_w1 = (factor_B - 1 + (int64_t)wt) / factor_B;
          int64_t max_w2 = factor_B * (int64_t)ws;
          int64_t min_w2 = (factor_A - 1 + (int64_t)ws) / factor_A;
#endif
          for (int o = m.l3; o < m.r3; o += 10) {
            int w1 = mark_answers[o + 0];
            int w2 = mark_answers[o + 1];
            T v = mark_answers[o + 2];
            T w = mark_answers[o + 3];
#ifdef WEIGHT
            if (w1 < min_w1 || w1 > max_w1 || w2 < min_w2 || w2 > max_w2) {
              continue;
            }
#endif
            // check duplication
            bool flag = false;
            for (int i = 1; i <= d1; i++) {
              if (st[i] == v || st[i] == w) {
                flag = true;
              }
            }
            if (flag) {
              continue;
            }
            ans_cnt++;
            memcpy(ans[d1 + d2], line, fixed_len[d1 + 1]);
            ans[d1 + d2] += str_len;
            memcpy(ans[d1 + d2], &mark_answers[o + 4], 24);
            ans[d1 + d2] += ((char *)&mark_answers[o + 4])[len_pos * 2];
            *(ans[d1 + d2] - 1) = '\n';
          }
        }
      }
    }
  }
}

template <typename T, bool mode>
inline void search_answers(ts_type ts[], Mark<T> marks[], T start, int pid,
                           int tid, int &ans_cnt, char *ans[]) {
  const Edge<T> *__restrict__ LA = (Edge<T> *)_LA;
  const size_t bs = maxl + 1;
  int ans_p = 0;
  register T st[5];
  char line[64];
  char *line_ptr = line;
  st[0] = start;
  memcpy(line_ptr, int_to_str[start], int_to_str_align);
  line_ptr += int_to_str[start][len_pos];
  if (mode) {
    const Edge2<T> *__restrict__ LA2 = (Edge2<T> *)_LA2;
    for (auto e = NA[start].L2; e < NA[start].R2; e++) {
      auto &v = LA2[e].im;
      auto &w = LA2[e].to;
      if (v <= start || w <= start) {
        continue;
      }
      st[1] = v;
      st[2] = w;
#ifdef WEIGHT
      int ws = LA2[e].w1;
      int wt = LA2[e].w2;
#else
      int ws = 0;
      int wt = 0;
#endif
      auto line_ptr_2 = line_ptr;
      memcpy(line_ptr_2, int_to_str[v], int_to_str_align);
      line_ptr_2 += int_to_str[v][len_pos];
      memcpy(line_ptr_2, int_to_str[w], int_to_str_align);
      line_ptr_2 += int_to_str[w][len_pos];
      if (unlikely(ts[w] == (ts_type)(start))) {
        print_answers<T, 2>(marks, st, start, w, pid, ans, ws, wt, ans_p, tid,
                            ans_cnt, line, line_ptr_2);
      }
#ifdef WEIGHT
      int64_t max_wu = factor_A * (int64_t)wt;
      int64_t min_wu = (factor_B - 1 + (int64_t)wt) / factor_B;
#endif
      for (auto f = NA[w].L2; f < NA[w].R2; f++) {
        auto &x = LA2[f].im;
        auto &y = LA2[f].to;
        // check later
        st[3] = x;
        st[4] = y;
#ifdef WEIGHT
        int wu = LA2[f].w1;
        int wv = LA2[f].w2;
#else
        int wu = 0;
        int wv = 0;
#endif
        if (unlikely(ts[y] == (ts_type)(start))) {
          if (x <= start || y <= start || x == v || y == v) {
            continue;
          }
#ifdef WEIGHT
          if (!check_weight(wt, wu)) {
            continue;
          }
#endif
          auto line_ptr_3 = line_ptr_2;
          memcpy(line_ptr_3, int_to_str[x], int_to_str_align);
          line_ptr_3 += int_to_str[x][len_pos];
          memcpy(line_ptr_3, int_to_str[y], int_to_str_align);
          line_ptr_3 += int_to_str[y][len_pos];
          print_answers<T, 4>(marks, st, start, y, pid, ans, ws, wv, ans_p, tid,
                              ans_cnt, line, line_ptr_3);
        }
      }
    }
  } else {
    for (auto e = NA[start].L; e < NA[start].R; e++) {
      auto &v = LA[e].to;
      if (v < start) {
        continue;
      }
      st[1] = v;
      for (auto f = NA[v].L; f < NA[v].R; f++) {
        auto &w = LA[f].to;
        if (w <= start) {
          continue;
        }
        st[2] = w;
#ifdef WEIGHT
        int ws = LA[e].w1;
        int wt = LA[f].w1;
        if (!check_weight(ws, wt)) {
          continue;
        }
#else
        int ws = 0;
        int wt = 0;
#endif
        auto line_ptr_2 = line_ptr;
        memcpy(line_ptr_2, int_to_str[v], int_to_str_align);
        line_ptr_2 += int_to_str[v][len_pos];
        memcpy(line_ptr_2, int_to_str[w], int_to_str_align);
        line_ptr_2 += int_to_str[w][len_pos];
        if (unlikely(ts[w] == (ts_type)(start))) {
          print_answers<T, 2>(marks, st, start, w, pid, ans, ws, wt, ans_p, tid,
                              ans_cnt, line, line_ptr_2);
        }
        for (auto g = NA[w].L; g < NA[w].R; g++) {
          auto &x = LA[g].to;
          st[3] = x;
          if (x < start) {
            continue;
          }
          for (auto h = NA[x].L; h < NA[x].R; h++) {
            auto &y = LA[h].to;
            st[4] = y;
            if (y <= start || y == w) {
              continue;
            }
#ifdef WEIGHT
            int wu = LA[g].w1;
            int wv = LA[h].w1;
            if (!check_weight(wu, wv)) {
              continue;
            }
#else
            int wu = 0;
            int wv = 0;
#endif
            if (unlikely(ts[y] == (ts_type)(start))) {
              if (x <= start || y <= start || x == v || y == v) {
                continue;
              }
#ifdef WEIGHT
              if (!check_weight(wt, wu)) {
                continue;
              }
#endif
              auto line_ptr_3 = line_ptr_2;
              memcpy(line_ptr_3, int_to_str[x], int_to_str_align);
              line_ptr_3 += int_to_str[x][len_pos];
              memcpy(line_ptr_3, int_to_str[y], int_to_str_align);
              line_ptr_3 += int_to_str[y][len_pos];
              print_answers<T, 4>(marks, st, start, y, pid, ans, ws, wv, ans_p,
                                  tid, ans_cnt, line, line_ptr_3);
            }
          }
        }
      }
    }
  }
}

template <typename T>
inline void run(int tid, T start, int pid, int &ans_cnt, char *ans[]) {
  Mark<T> *marks = (Mark<T> *)_marks[tid];
  ts_type *ts = _ts[tid];
  if (unlikely(start == 0)) {
    for (int i = 0; i < nodes_count; i++) {
      ts[i] = -1;
    }
    for (int i = 0; i < nodes_count; i++) {
      marks[i].s1 = marks[i].s2 = marks[i].s3 = -1;
    }
  }
  Timer::start(SP_STEP1);
  if (adj_list_2_mode & 2) {
    process_marks<T, true>(ts, marks, start);
  } else {
    process_marks<T, false>(ts, marks, start);
  }
  Timer::end(SP_STEP1);
  Timer::start(SP_STEP2);
  if (adj_list_2_mode & 1) {
    search_answers<T, true>(ts, marks, start, pid, tid, ans_cnt, ans);
  } else {
    search_answers<T, false>(ts, marks, start, pid, tid, ans_cnt, ans);
  }
  Timer::end(SP_STEP2);
}

template <typename T>
void worker(int tid) {
  int ans_cnt = 0;
  int pid;  // part_id
  const int bs = maxl + 1;
  for (int i = 0; i < bs; i++) {
    answers_ptr_buf[tid][i] = answers_initial_pointer(tid, i);
  }
  while ((pid = get_job()) < num_parts) {
    // fprintf(stderr, "worker %d running part %d\n", tid, pid);
    int start = (pid + 0) * nodes_count / num_parts;
    int end = (pid + 1) * nodes_count / num_parts;
    if (start == end) {
      continue;
    }
    register char *ans[bs];
    for (int i = 0; i < bs; i++) {
      ans[i] = answers_ptr[pid][i].first = answers_ptr_buf[tid][i];
    }
    for (int i = start; i < end; i++) {
      run<T>(tid, i, pid, ans_cnt, ans);
      finished[i] = true;
    }
    for (int i = 0; i < bs; i++) {
      answers_ptr[pid][i].second = answers_ptr_buf[tid][i] = ans[i];
    }
    update_graph<T>();
  }
  std::atomic_fetch_add_explicit(&answers_count, ans_cnt,
                                 std::memory_order_relaxed);
}

namespace TPP {

static std::atomic<int> deg[maxn][2];
static int pt[maxn][2][2];  // {u}{[ST]}{[LR]}
static int id_new[maxn];

inline void topo_prune_recursive(int tid) {
  static int Q[num_threads][maxn];
  for (int s = tid; s < nodes_count; s += num_threads) {
    for (int dir = 0; dir <= 1; dir++) {
      if (deg[s][dir] == 0) {
        int L = 0, R = 1;
        deg[s][dir] = -1;
        Q[tid][0] = s;
        while (L < R) {
          int u = Q[tid][L++];
          for (int i = pt[u][dir ^ 1][0]; i < pt[u][dir ^ 1][1]; i++) {
            auto v = *(&edges_list[0][0] + (ii_list[i].ofs ^ 1));
            if (std::atomic_fetch_sub_explicit(
                    &deg[v][dir], 1, std::memory_order_relaxed) == 1) {
              deg[v][dir] = -1;
              Q[tid][R++] = v;
            }
          }
        }
      }
    }
  }
}

inline void calculate_deg(int tid) {
  int start = tid * edges_count * 2 / num_threads;
  int end = (tid + 1) * edges_count * 2 / num_threads;
  if (start == end) {  // must check
    return;
  }
  while (start > 0 && start < edges_count * 2 &&
         (ii_list[start - 1].val) == (ii_list[start].val)) {
    start++;
  }
  while (end > 0 && end < edges_count * 2 &&
         (ii_list[end - 1].val) == (ii_list[end].val)) {
    end++;
  }
  for (int i = start; i < end; i++) {
    const auto ii_val = ii_list[i].val;
    const auto val = ii_val >> 1;
    const auto dir = ii_val & 1;
    if (pt[val][dir][0] == 0) {
      pt[val][dir][0] = i;
    }
    pt[val][dir][1] = i + 1;
    deg[val][dir]++;
  }
}

inline void relabel() {
  int cnt = 0;
  for (int i = 0; i < nodes_count; i++) {
    if (deg[i][0] > 0 && deg[i][1] > 0) {
      id_rev[cnt] = id_rev[i];
      id_new[i] = cnt++;
    } else {
      id_new[i] = -1;
    }
  }
  nodes_count = cnt;
  int tot = 0;
  for (int i = 0; i < edges_count; i++) {
    int u = edges_list[i][0];
    int v = edges_list[i][1];
    int w = edges_list[i][2];
    if (id_new[u] != -1 && id_new[v] != -1) {
      edges_list[tot][0] = id_new[u];
      edges_list[tot][1] = id_new[v];
      edges_list[tot][2] = w;
      tot++;
    }
  }
  edges_count = tot;
}

inline void topo_prune() {
  Timer::start(SP_DEG);
  if (topo_prune_mode != 0) {
    parallel_run<num_threads>(calculate_deg);
  }
  Timer::end(SP_DEG);
  Timer::start(SP_REM);
  if (topo_prune_mode == 2) {
    parallel_run<num_threads>(topo_prune_recursive);
  }
  Timer::end(SP_REM);
  Timer::start(SP_RLB);
  if (topo_prune_mode != 0) {
    relabel();
  }
  Timer::end(SP_RLB);
}
};  // namespace TPP

inline void ID_LSH() {  // need parallelization
  Timer::start(SP_PSRS);
  PSRS::parasort(edges_count * 2, ii_list);
  Timer::end(SP_PSRS);
  Timer::start(SP_UNIQ);
  nodes_count = 0;
  int j = 0;
  int tmp = 0;
  for (int i = 0; i < edges_count * 2; i++) {
    const auto &val = ii_list[i].val >> 1;
    const auto &ofs = ii_list[i].ofs;
    if (i == 0 || tmp < val) {
      id_rev[j++] = val;
    }
    tmp = val;
    *(&edges_list[0][0] + ofs) = j - 1;
    ii_list[i].val = ((j - 1) << 1) | (ofs & 1);
  }
  nodes_count = j;
  Timer::end(SP_UNIQ);
}

template <typename T>
void prepare_DS_0(int tid) {  // int to str
  int start = tid * nodes_count / num_threads;
  int end = (tid + 1) * nodes_count / num_threads;
  for (int i = start; i < end; i++) {
#ifdef USE_LSH
    int x = id_rev[i];
#else
    int x = i;
#endif
    int p = 0;
    if (x == 0) {
      int_to_str[i][p++] = '0';
    } else {
      while (x) {
        int_to_str[i][p++] = x % 10 + '0';
        x /= 10;
      }
    }
    std::reverse(int_to_str[i], int_to_str[i] + p);
    int_to_str[i][p++] = cm;
    int_to_str[i][len_pos] = p;
  }
}

template <typename T>
void prepare_DS_1(int tid) {
  for (int i = 0; i < edges_count; i++) {
    T u = edges_list[i][0];
    T v = edges_list[i][1];
    if (my_hash<T, num_threads>(u) == tid) {
      NA[u].C++;
    }
    if (my_hash<T, num_threads>(v) == tid) {
      NB[v].C++;
    }
  }
}

template <typename T>
void prepare_DS_2(int tid) {
  for (int i = 0; i < edges_count; i++) {
    T u = edges_list[i][0];
    T v = edges_list[i][1];
    if (my_hash<T, num_threads>(u) == tid) {
      NA[u].C2 += NA[v].C;
    }
    if (my_hash<T, num_threads>(v) == tid) {
      NB[v].C2 += NB[u].C;
    }
  }
}

template <typename T>
void prepare_DS_3(int tid) {
  if (num_threads == 1 || tid == 0) {
    int P = 0;
    int P2 = 0;
    for (int i = 0; i < nodes_count; i++) {
      NA[i].R = NA[i].L = P;
      NA[i].R2 = NA[i].L2 = P2;
      P += NA[i].C;
      P2 += NA[i].C2;
    }
  }
  if (num_threads == 1 || tid == 1) {
    int P = 0;
    int P2 = 0;
    for (int i = 0; i < nodes_count; i++) {
      NB[i].R = NB[i].L = P;
      NB[i].R2 = NB[i].L2 = P2;
      P += NB[i].C;
      P2 += NB[i].C2;
    }
  }
}

template <typename T>
void prepare_DS_4(int tid) {
  Edge<T> *__restrict__ LA = (Edge<T> *)_LA;
  Edge<T> *__restrict__ LB = (Edge<T> *)_LB;
  for (int i = 0; i < edges_count; i++) {
    T u = edges_list[i][0];
    T v = edges_list[i][1];
    int w = edges_list[i][2];
    if (my_hash<T, num_threads>(u) == tid) {
#ifdef WEIGHT
      LA[NA[u].R++] = {v, w};
#else
      LA[NA[u].R++] = {v};
#endif
    }
    if (my_hash<T, num_threads>(v) == tid) {
#ifdef WEIGHT
      LB[NB[v].R++] = {u, w};
#else
      LB[NB[v].R++] = {u};
#endif
    }
  }
}

template <typename T>
void prepare_DS_5(int tid) {
  Edge<T> *__restrict__ LA = (Edge<T> *)_LA;
  Edge<T> *__restrict__ LB = (Edge<T> *)_LB;
  for (int i = 0; i < nodes_count; i++) {
    if (my_hash<T, num_threads>(i) == tid) {
      std::sort(&LA[NA[i].L], &LA[NA[i].R]);
      std::sort(&LB[NB[i].L], &LB[NB[i].R]);
    }
  }
}

template <typename T>
void prepare_DS_6(int tid) {
  static std::atomic<int> current_job;
  Edge<T> *__restrict__ LA = (Edge<T> *)_LA;
  Edge2<T> *__restrict__ LA2 = (Edge2<T> *)_LA2;
  int pid;
  if (!(adj_list_2_mode & 1)) {
    return;
  }
  while ((pid = std::atomic_fetch_add_explicit(
              &current_job, 1, std::memory_order_relaxed)) < num_parts) {
    int start = (pid + 0) * nodes_count / num_parts;
    int end = (pid + 1) * nodes_count / num_parts;
    for (int u = start; u < end; u++) {
      for (auto e = NA[u].L; e < NA[u].R; e++) {
        auto &v = LA[e].to;
        auto &w1 = LA[e].w1;
#ifdef WEIGHT
        int64_t max_w2 = factor_A * (int64_t)w1;
        int64_t min_w2 = (factor_B - 1 + (int64_t)w1) / factor_B;
#endif
        for (auto f = NA[v].L; f < NA[v].R; f++) {
          auto &w = LA[f].to;
          if (w != u) {
#ifdef WEIGHT
            auto &w2 = LA[f].w1;
            if (w2 < min_w2 || w2 > max_w2) {
              continue;
            }
#endif
            auto &z = LA2[NA[u].R2++];
            z.to = w;
            z.im = v;
#ifdef WEIGHT
            z.w1 = w1;
            z.w2 = w2;
#endif
          }
        }
      }
    }
  }
}

template <typename T>
void prepare_DS_7(int tid) {
  static std::atomic<int> current_job;
  Edge<T> *__restrict__ LB = (Edge<T> *)_LB;
  Edge2<T> *__restrict__ LB2 = (Edge2<T> *)_LB2;
  int pid;
  if (!(adj_list_2_mode & 2)) {
    return;
  }
  while ((pid = std::atomic_fetch_add_explicit(
              &current_job, 1, std::memory_order_relaxed)) < num_parts) {
    int start = (pid + 0) * nodes_count / num_parts;
    int end = (pid + 1) * nodes_count / num_parts;
    for (int u = start; u < end; u++) {
      if (adj_list_2_mode & 2) {
        for (auto e = NB[u].L; e < NB[u].R; e++) {
          auto &v = LB[e].to;
          auto &w1 = LB[e].w1;
#ifdef WEIGHT
          int64_t max_w2 = factor_B * (int64_t)w1;
          int64_t min_w2 = (factor_A - 1 + (int64_t)w1) / factor_A;
#endif
          for (auto f = NB[v].L; f < NB[v].R; f++) {
            auto &w = LB[f].to;
            if (w != u) {
#ifdef WEIGHT
              auto &w2 = LB[f].w1;
              if (w2 < min_w2 || w2 > max_w2) {
                continue;
              }
#endif
              auto &z = LB2[NB[u].R2++];
              z.to = w;
              z.im = v;
#ifdef WEIGHT
              z.w1 = w1;
              z.w2 = w2;
#endif
            }
          }
        }
      }
    }
  }
}

template <typename T>
void prepare_DS() {
  {
    Timer::start(SP_DS_0);
    parallel_run<num_threads>(prepare_DS_0<T>);
    Timer::end(SP_DS_0);
  }
  {
    Timer::start(SP_DS_1);
    parallel_run<num_threads>(prepare_DS_1<T>);
    Timer::end(SP_DS_1);
  }
  {
    Timer::start(SP_DS_2);
    parallel_run<num_threads>(prepare_DS_2<T>);
    Timer::end(SP_DS_2);
  }
  {
    Timer::start(SP_DS_3);
    parallel_run<num_threads>(prepare_DS_3<T>);
    Timer::end(SP_DS_3);
  }
  {
    Timer::start(SP_DS_4);
    parallel_run<num_threads>(prepare_DS_4<T>);
    Timer::end(SP_DS_4);
  }
  {
    Timer::start(SP_DS_5);
    parallel_run<num_threads>(prepare_DS_5<T>);
    Timer::end(SP_DS_5);
  }
  {
    Timer::start(SP_DS_6);
    parallel_run<num_threads>(prepare_DS_6<T>);
    Timer::end(SP_DS_6);
  }
  {
    Timer::start(SP_DS_7);
    parallel_run<num_threads>(prepare_DS_7<T>);
    Timer::end(SP_DS_7);
  }
}

template <int num_tasks>
void reader(int tid) {
  static std::atomic<int> current_job;
  int start = tid * input_file_size / num_tasks;
  int end = (tid + 1) * input_file_size / num_tasks;
  while (start > 0 && start < input_file_size && region[start] != '\n') {
    start++;
  }
  while (end > 0 && end < input_file_size && region[end] != '\n') {
    end++;
  }
  if (start == end) {
    return;
  }
  int o = start;
  char c;
  while (o < end && ((c = region[o]), !is_digit(c))) {
    o++;
  }
  int basket_id = std::atomic_fetch_add_explicit(&current_job, 1,
                                                 std::memory_order_relaxed);
  int cnt = 0;
  while (o < end) {
    int flag = 0;
    int a[3] = {};
    for (int i = 0; i < 3; i++) {
      while ((c = region[o++]), is_digit(c)) {
        a[i] *= 10;
        a[i] += c - '0';
      }  // assume seperated by a comma
    }
    while (o < end && ((c = region[o]), !is_digit(c))) {
      o++;
    }  // read LF or CRLF
    if (a[2] == 0) {
      continue;
    }
    int idx = basket_id * basket_size + cnt;
    cnt++;
    if (cnt == basket_size) {
      cnt = 0;
      basket_id = std::atomic_fetch_add_explicit(&current_job, 1,
                                                 std::memory_order_relaxed);
    }
    edges_list[idx][0] = a[0];
    edges_list[idx][1] = a[1];
    edges_list[idx][2] = a[2];
#ifdef USE_LSH
    ii_list[idx * 2 + 0] = {(uint32_t)(&edges_list[idx][0] - &edges_list[0][0]),
                            (((uint32_t)a[0]) << 1) | 0};
    ii_list[idx * 2 + 1] = {(uint32_t)(&edges_list[idx][1] - &edges_list[0][0]),
                            (((uint32_t)a[1]) << 1) | 1};
#else
    nodes_count = std::max(nodes_count, std::max(a[0], a[1]) + 1);
#endif
  }
  edges_count_each[tid] = basket_id * basket_size + cnt;
}

void input() {
  int u, v, w;
#ifdef USE_MMAP_INPUT
  size_t map_len;
#endif
  {
    FILE *pFILE = fopen(input_file_name, "r");
    fseek(pFILE, 0, SEEK_END);
    auto file_size = input_file_size = ftell(pFILE);
    fclose(pFILE);
    in_fd = open(input_file_name, O_RDONLY);
#ifdef USE_MMAP_INPUT
    {
      map_len = file_size;
      region = (char *)(mmap(NULL, map_len, PROT_READ,
                             MAP_FILE | MAP_PRIVATE | MAP_NORESERVE, in_fd, 0));
      madvise(region, map_len, MADV_WILLNEED);
    }
#else
    {
      region = (char *)malloc(file_size);
      read(in_fd, region, file_size);
    }
#endif
    const int num_tasks = num_threads;
    edges_count_each.resize(num_tasks);
    parallel_run<num_tasks>(reader<num_tasks>);
    std::sort(edges_count_each.begin(), edges_count_each.end());
    edges_count = edges_count_each.back();
    for (int i = 0; i < num_tasks; i++) {
      int start = edges_count_each[i];
      do {
        while (edges_count > 1 && edges_list[edges_count - 1][0] ==
                                      edges_list[edges_count - 1][1]) {
          edges_count--;
        }
        if (start >= edges_count - 1) {
          break;
        }
        edges_count--;
#ifdef USE_LSH
        ii_list[start * 2 + 0] = {
            (uint32_t)(&edges_list[start][0] - &edges_list[0][0]),
            (((uint32_t)edges_list[edges_count][0]) << 1) | 0};
        ii_list[start * 2 + 1] = {
            (uint32_t)(&edges_list[start][1] - &edges_list[0][0]),
            (((uint32_t)edges_list[edges_count][1]) << 1) | 1};
#endif
        memcpy(&edges_list[start][0], &edges_list[edges_count][0],
               sizeof(edges_list[0][0]) * 4);
        start++;
      } while (start % basket_size != 0);
    }
#ifdef USE_MMAP_INPUT
    { munmap(region, map_len); }
#else
    { free(region); }
#endif
    close(in_fd);
  }
}

template <typename T>
void output() {
  const int bs = maxl + 1;
  const size_t ps = getpagesize();
  const int num_tasks = num_threads;
  int tot = answers_count.load();
  size_t text_size = 0;
  char buf[16] = {};
  sprintf(buf, "%d\n", tot);
  size_t buf_len = strlen(buf);
  text_size += buf_len;
  for (int tid = 0; tid < num_threads; tid++) {
    for (int len = minl; len <= maxl; len++) {
      text_size +=
          answers_ptr_buf[tid][len] - answers_initial_pointer(tid, len);
    }
  }
#ifdef INFO
  Color::green();
  fprintf(stderr, "NUM_ANSWERS = %d\n", tot);
  fprintf(stderr, "TEXT_SIZE = %zu bytes\n", text_size);
  Color::reset();
#endif
  out_fd = open(output_file_name, O_RDWR | O_CREAT | O_NONBLOCK, (mode_t)0666);
  ftruncate(out_fd, text_size);
  pwrite(out_fd, buf, buf_len, 0);
  size_t offset = buf_len;
  for (int len = minl; len <= maxl; len++) {
    for (int i = 0; i < num_parts; i++) {
      size_t sz = answers_ptr[i][len].second - answers_ptr[i][len].first;
      pwrite(out_fd, answers_ptr[i][len].first, sz, offset);
      offset += sz;
    }
  }
  close(out_fd);
}

void print_info() {
  Color::yellow();
  fprintf(stderr, "---------------------------------\n");
  Color::reset();
  Color::blue();
  Timer::print(SP_INPUT, "INPUT");
#ifdef USE_LSH
  Color::cyan();
  Timer::print(SP_LSH, "LSH");
  Color::green();
  Timer::print(SP_PSRS, "PSRS");
  Timer::print(SP_UNIQ, "UNIQ");
  Color::reset();
#endif
#ifdef ALLOW_TOPO_PRUNE
  Color::cyan();
  Timer::print(SP_PRUNE, "PRUNE");
  Color::green();
  Timer::print(SP_DEG, "DEG");
  Timer::print(SP_REM, "REM");
  Timer::print(SP_RLB, "RLB");
  Color::reset();
#endif
  Color::cyan();
  Timer::print(SP_DS, "DS");
  Color::green();
  Timer::print(SP_DS_0, "DS_0");
  Timer::print(SP_DS_1, "DS_1");
  Timer::print(SP_DS_2, "DS_2");
  Timer::print(SP_DS_3, "DS_3");
  Timer::print(SP_DS_4, "DS_4");
  Timer::print(SP_DS_5, "DS_5");
  Timer::print(SP_DS_6, "DS_6");
  Timer::print(SP_DS_7, "DS_7");
  Color::magenta();
  Timer::print(SP_RUN, "RUN");
  Color::green();
  Timer::print(SP_STEP1, "STEP1");
  Timer::print(SP_STEP2, "STEP2");
  Timer::print(SP_UPD, "UPD");
  Color::blue();
  Timer::print(SP_OUTPUT, "OUTPUT");
  Color::orange();
  Timer::print(SP_TOTAL, "TOTAL");
  Color::reset();
  Color::yellow();
  fprintf(stderr, "---------------------------------\n");
  Color::reset();
#ifdef STATISTICS
  Color::reset();
  fprintf(stderr, "SEARCH_COUNT = %lld\n", search_count.load());
  fprintf(stderr, "HIT_COUNT = %lld\n", hit_count.load());
  fprintf(stderr, "FB_COUNT = %lld\n", fb_count.load());
  fprintf(stderr, "SEARCH_EFFICIENCY = %.6lf\n",
          (double)hit_count.load() / search_count.load());
  fprintf(stderr, "FB_RATIO = %.6lf\n",
          (double)fb_count.load() / hit_count.load());
  Color::reset();
#endif
}

int main() {
#ifdef KASHI
  auto start_time = std::chrono::clock_type::now();
#endif
#ifdef ALLOW_UINT16
  bool use_uint16_t = true;
#endif
  if (true) {
    delay<pre_sleep_time>();
  }
  Timer::start(SP_TOTAL);
  Timer::start(SP_INPUT);
  input();
  Timer::end(SP_INPUT);
#ifdef USE_LSH
  Timer::start(SP_LSH);
  ID_LSH();
  Timer::end(SP_LSH);
#endif
#ifdef ALLOW_TOPO_PRUNE
  Timer::start(SP_PRUNE);
  TPP::topo_prune();
  Timer::end(SP_PRUNE);
#endif
#ifdef ALLOW_UINT16
  if (nodes_count >= UINT16_MAX) {  // choose uint16_t or uint32_t
    use_uint16_t = false;
  } else {
    use_uint16_t = true;
  }
#endif
  Timer::start(SP_DS);
#ifdef ALLOW_UINT16
  if (use_uint16_t) {
    prepare_DS<uint16_t>();
  } else {
    prepare_DS<uint32_t>();
  }
#else
  { prepare_DS<uint32_t>(); }
#endif
  Timer::end(SP_DS);
#ifdef INFO
  Color::reset();
  Color::yellow();
  fprintf(stderr, "topological_prune_mode = %d\n", topo_prune_mode);
  fprintf(stderr, "adjacent_list_2_mode = %d\n", adj_list_2_mode);
  Color::green();
  fprintf(stderr, "graph_nodes_count = %d\n", (int)nodes_count);
  fprintf(stderr, "graph_edges_count = %d\n", (int)edges_count);
  fprintf(stderr, "forward_star_table_2_size (orig) = %d\n",
          NA[nodes_count - 1].R2);
  fprintf(stderr, "forward_star_table_2_size (rev)  = %d\n",
          NB[nodes_count - 1].R2);
  Color::reset();
#endif
#ifdef INFO
  Color::yellow();
#ifdef ALLOW_UINT16
  fprintf(stderr, "Using %s\n", use_uint16_t ? "uint16_t" : "uint32_t");
#else
  fprintf(stderr, "Using %s\n", "uint32_t");
#endif
#ifdef USE_MMAP
  fprintf(stderr, "Using mmap()\n");
#else
  fprintf(stderr, "Using write()\n");
#endif
  Color::reset();
  Color::green();
  fprintf(stderr, "NUM_THREADS = %d\n", num_threads);
  fprintf(stderr, "NUM_PARTS = %d\n", num_parts);
  Color::reset();
#endif
  Timer::start(SP_RUN);
#ifdef ALLOW_UINT16
  if (use_uint16_t) {
    init_answers<uint16_t>();
  } else {
    init_answers<uint32_t>();
  }
#else
  { init_answers<uint32_t>(); }
#endif
#ifdef ALLOW_UINT16
  if (use_uint16_t) {
    parallel_run<num_threads>(worker<uint16_t>);
  } else {
    parallel_run<num_threads>(worker<uint32_t>);
  }
#else
  { parallel_run<num_threads>(worker<uint32_t>); }
#endif
  Timer::end(SP_RUN);
  Timer::start(SP_OUTPUT);
#ifdef ALLOW_UINT16
  if (use_uint16_t) {
    output<uint16_t>();
  } else {
    output<uint32_t>();
  }
#else
  { output<uint32_t>(); }
#endif
  Timer::end(SP_OUTPUT);
  Timer::end(SP_TOTAL);
#ifdef INFO
  print_info();
#endif
  delay<post_sleep_time>();
#ifdef KASHI
  const uint64_t target_time = Timer::gettime(SP_OUTPUT) * 100;
  while (true) {
    auto end_time = std::chrono::clock_type::now();
    auto duration = end_time - start_time;
    uint64_t secs =
        std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count();
    if (secs >= target_time) {
      break;
    }
  }
#endif
  // exit(EXIT_SUCCESS);
  _Exit(EXIT_SUCCESS);
}