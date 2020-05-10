#include <bits/stdc++.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>

using namespace std;

inline const char* GetDigitsLut() {
  static const char cDigitsLut[200] = {
      '0', '0', '0', '1', '0', '2', '0', '3', '0', '4', '0', '5', '0', '6', '0',
      '7', '0', '8', '0', '9', '1', '0', '1', '1', '1', '2', '1', '3', '1', '4',
      '1', '5', '1', '6', '1', '7', '1', '8', '1', '9', '2', '0', '2', '1', '2',
      '2', '2', '3', '2', '4', '2', '5', '2', '6', '2', '7', '2', '8', '2', '9',
      '3', '0', '3', '1', '3', '2', '3', '3', '3', '4', '3', '5', '3', '6', '3',
      '7', '3', '8', '3', '9', '4', '0', '4', '1', '4', '2', '4', '3', '4', '4',
      '4', '5', '4', '6', '4', '7', '4', '8', '4', '9', '5', '0', '5', '1', '5',
      '2', '5', '3', '5', '4', '5', '5', '5', '6', '5', '7', '5', '8', '5', '9',
      '6', '0', '6', '1', '6', '2', '6', '3', '6', '4', '6', '5', '6', '6', '6',
      '7', '6', '8', '6', '9', '7', '0', '7', '1', '7', '2', '7', '3', '7', '4',
      '7', '5', '7', '6', '7', '7', '7', '8', '7', '9', '8', '0', '8', '1', '8',
      '2', '8', '3', '8', '4', '8', '5', '8', '6', '8', '7', '8', '8', '8', '9',
      '9', '0', '9', '1', '9', '2', '9', '3', '9', '4', '9', '5', '9', '6', '9',
      '7', '9', '8', '9', '9'};
  return cDigitsLut;
}

inline int int2str(int value, char*& buffer) {
  const char* cDigitsLut = GetDigitsLut();
  static int sz = 0;

  if (value < 10)
    sz = 1;
  else if (value < 100)
    sz = 2;
  else if (value < 1000)
    sz = 3;
  else if (value < 10000)
    sz = 4;
  else if (value < 100000)
    sz = 5;
  else if (value < 1000000)
    sz = 6;
  else if (value < 10000000)
    sz = 7;
  else if (value < 100000000)
    sz = 8;
  else if (value < 1000000000)
    sz = 9;
  else
    sz = 10;

  if (value < 10000) {
    const uint32_t d1 = (value / 100) << 1;
    const uint32_t d2 = (value % 100) << 1;

    if (value >= 1000) *buffer++ = cDigitsLut[d1];
    if (value >= 100) *buffer++ = cDigitsLut[d1 + 1];
    if (value >= 10) *buffer++ = cDigitsLut[d2];
    *buffer++ = cDigitsLut[d2 + 1];
  } else if (value < 100000000) {
    // value = bbbbcccc
    const uint32_t b = value / 10000;
    const uint32_t c = value % 10000;

    const uint32_t d1 = (b / 100) << 1;
    const uint32_t d2 = (b % 100) << 1;

    const uint32_t d3 = (c / 100) << 1;
    const uint32_t d4 = (c % 100) << 1;

    if (value >= 10000000) *buffer++ = cDigitsLut[d1];
    if (value >= 1000000) *buffer++ = cDigitsLut[d1 + 1];
    if (value >= 100000) *buffer++ = cDigitsLut[d2];
    *buffer++ = cDigitsLut[d2 + 1];

    *buffer++ = cDigitsLut[d3];
    *buffer++ = cDigitsLut[d3 + 1];
    *buffer++ = cDigitsLut[d4];
    *buffer++ = cDigitsLut[d4 + 1];
  } else {
    // value = aabbbbcccc in decimal

    const uint32_t a = value / 100000000;  // 1 to 42
    value %= 100000000;

    if (a >= 10) {
      const unsigned i = a << 1;
      *buffer++ = cDigitsLut[i];
      *buffer++ = cDigitsLut[i + 1];
    } else
      *buffer++ = static_cast<char>('0' + static_cast<char>(a));

    const uint32_t b = value / 10000;  // 0 to 9999
    const uint32_t c = value % 10000;  // 0 to 9999

    const uint32_t d1 = (b / 100) << 1;
    const uint32_t d2 = (b % 100) << 1;

    const uint32_t d3 = (c / 100) << 1;
    const uint32_t d4 = (c % 100) << 1;

    *buffer++ = cDigitsLut[d1];
    *buffer++ = cDigitsLut[d1 + 1];
    *buffer++ = cDigitsLut[d2];
    *buffer++ = cDigitsLut[d2 + 1];
    *buffer++ = cDigitsLut[d3];
    *buffer++ = cDigitsLut[d3 + 1];
    *buffer++ = cDigitsLut[d4];
    *buffer++ = cDigitsLut[d4 + 1];
  }

  *buffer = ',';
  return sz + 1;
}

// inline static int int2str(int x, char* &it) {
//    static const char digits[] =
//            "0001020304050607080910111213141516171819"
//            "2021222324252627282930313233343536373839"
//            "4041424344454647484950515253545556575859"
//            "6061626364656667686970717273747576777879"
//            "8081828384858687888990919293949596979899";
//    static int __sz, __nxt;
//    if (x < 10) __sz = 1;
//    else if (x < 100) __sz = 2;
//    else if (x < 1000) __sz = 3;
//    else if (x < 10000) __sz = 4;
//    else if (x < 100000) __sz = 5;
//    else if (x < 1000000) __sz = 6;
//    else __sz = 7;
//    __nxt = __sz - 1;
//
//    while(x >= 100) {
//        const int i = (x % 100) << 1;
//        x /= 100;
//        it[__nxt - 1] = digits[i];
//        it[__nxt] = digits[i + 1];
//        __nxt -= 2;
//    }
//
//    if(x < 10) it[__nxt] = x ^ 48;
//    else {
//        const int i = x << 1;
//        it[__nxt - 1] = digits[i];
//        it[__nxt] = digits[i + 1];
//    }
//
//    it[__sz] = ',';
//    return __sz + 1;
//}

uint8_t MAXD = 100;
#define MAX_RATIO 3.0
#define MIN_RATIO 0.2

#define ID_MAX_NUM 4000000

struct Vertex {
  int v;
  int w;
  Vertex() {}
  Vertex(int& v, int& w) : v(v), w(w) {}
  inline bool operator<(const Vertex& ver) const { return v > ver.v; }
};

// int * d_c;
unordered_map<int, int> d_c;
// map<int, map<int, float>> weights;

int n, N, all_nodes[2 * ID_MAX_NUM], edges[3 * ID_MAX_NUM];
Vertex* AK;
Vertex* AK_reverse;
uint8_t in_deg[ID_MAX_NUM], out_deg[ID_MAX_NUM];
uint8_t AK_cnt[ID_MAX_NUM], AK_reverse_cnt[ID_MAX_NUM];

int* res3_t1 = new int[20000000 * 3];
int* res4_t1 = new int[20000000 * 4];
int* res5_t1 = new int[20000000 * 5];
int* res6_t1 = new int[20000000 * 6];
int* res7_t1 = new int[20000000 * 7];

int* res3_t2 = new int[20000000 * 3];
int* res4_t2 = new int[20000000 * 4];
int* res5_t2 = new int[20000000 * 5];
int* res6_t2 = new int[20000000 * 6];
int* res7_t2 = new int[20000000 * 7];

int* res3_t3 = new int[20000000 * 3];
int* res4_t3 = new int[20000000 * 4];
int* res5_t3 = new int[20000000 * 5];
int* res6_t3 = new int[20000000 * 6];
int* res7_t3 = new int[20000000 * 7];

int* res3_t4 = new int[20000000 * 3];
int* res4_t4 = new int[20000000 * 4];
int* res5_t4 = new int[20000000 * 5];
int* res6_t4 = new int[20000000 * 6];
int* res7_t4 = new int[20000000 * 7];

int res_cnt[4][5];

char* answer = new char[7 * 10 * 20000000];

char str[ID_MAX_NUM * 11];
char* map_int_str[ID_MAX_NUM];
uint8_t map_int_len[ID_MAX_NUM];

int* res3[4] = {res3_t1, res3_t2, res3_t3, res3_t4};
int* res4[4] = {res4_t1, res4_t2, res4_t3, res4_t4};
int* res5[4] = {res5_t1, res5_t2, res5_t3, res5_t4};
int* res6[4] = {res6_t1, res6_t2, res6_t3, res6_t4};
int* res7[4] = {res7_t1, res7_t2, res7_t3, res7_t4};

void TopoSortRemove() {
  queue<int> q;

  // remove in
  for (int i = 0; i < n; ++i) {
    if (!in_deg[i]) q.push(i);
  }

  while (!q.empty()) {
    int u = q.front();
    q.pop();
    Vertex* g = AK + u * MAXD;
    for (int i = 0; i < AK_cnt[u]; ++i) {
      if (--in_deg[g->v] == 0) q.push(g->v);
      ++g;
    }
  }

  for (int i = 0; i < n; ++i) {
    if (!in_deg[i]) {
      AK_cnt[i] = 0;
      AK_reverse_cnt[i] = 0;
    }
  }

  // remove out
  for (int i = 0; i < n; ++i) {
    if (!out_deg[i]) q.push(i);
  }

  while (!q.empty()) {
    int u = q.front();
    q.pop();
    Vertex* r = AK_reverse + u * MAXD;
    for (int i = 0; i < AK_reverse_cnt[u]; ++i) {
      if (--out_deg[r->v] == 0) q.push(r->v);
      ++r;
    }
  }

  for (int i = 0; i < n; ++i) {
    if (!out_deg[i]) {
      AK_cnt[i] = 0;
      AK_reverse_cnt[i] = 0;
    }
  }

  while (AK_cnt[n - 1] == 0) --n;
}

int CountSort(int* vec, int n, int max_val) {
  int* ele_cnt = new int[max_val + 1];
  int i = 0;
  while (i < n) {
    ele_cnt[vec[i]] = 1;
    ++i;
  }
  i = 0;
  for (int j = 0; j <= max_val; j++) {
    if (ele_cnt[j]) {
      vec[i++] = j;
    }
  }
  return i;
}

void MyLoadData(const char* inputFile) {
#ifdef TEST
  auto start_time = chrono::steady_clock::now();
#endif

  int fd = open(inputFile, O_RDONLY);

  struct stat info;
  stat(inputFile, &info);
  size_t buf_len = info.st_size;

  //    cout << buf_len << endl;
  //    int buf_len = 10 * 1024 * 1204;
  char* buf = (char*)mmap(NULL, buf_len, PROT_READ, MAP_PRIVATE, fd, 0);

  int start, end, weight;
  N = 0;

  int edge_cnt = 0;

  do {
    start = 0;
    while (*buf != ',') {
      start = (start << 3) + (start << 1) + (*buf ^ 48);
      buf++;
    }

    buf++;

    end = 0;
    while (*buf != ',') {
      end = (end << 3) + (end << 1) + (*buf ^ 48);
      buf++;
    }

    buf++;

    //        cout << start << "," << end << "," << N << endl;

    weight = 0;
    while (*buf != '\r' && *buf != '\n') {
      weight = (weight << 3) + (weight << 1) + (*buf ^ 48);
      buf++;
    }

    if (*buf == '\r') buf++;
    buf++;

    all_nodes[n++] = start;
    all_nodes[n++] = end;

    edges[edge_cnt++] = start;
    edges[edge_cnt++] = end;
    edges[edge_cnt++] = weight;

    N = max(N, start);
    N = max(N, end);

  } while (*buf != '\0');

  //    sort(all_nodes, all_nodes + n);
  sort(all_nodes, all_nodes + n);
  n = unique(all_nodes, all_nodes + n) - all_nodes;
  //    n = CountSort(all_nodes, n, N);

  //    d_c = new int[N + 1];
  char* it = str;
  for (int i = 0; i < n; i++) {
    d_c[all_nodes[i]] = i;
    map_int_str[i] = it;
    map_int_len[i] = int2str(all_nodes[i], it);
    it += 11;
  }

  for (int i = 0; i < edge_cnt; i += 3) {
    start = d_c[edges[i]];
    end = d_c[edges[i + 1]];
    ++out_deg[start];
    ++in_deg[end];
    MAXD = max(MAXD, out_deg[start]);
    MAXD = max(MAXD, in_deg[end]);
  }

  AK = new Vertex[ID_MAX_NUM * MAXD];
  AK_reverse = new Vertex[ID_MAX_NUM * MAXD];

  Vertex* pv;
  for (int i = 0; i < edge_cnt; i += 3) {
    start = d_c[edges[i]];
    end = d_c[edges[i + 1]];
    weight = edges[i + 2];

    pv = AK + start * MAXD + AK_cnt[start];
    pv->v = end;
    pv->w = weight;
    AK_cnt[start]++;

    pv = AK_reverse + end * MAXD + AK_reverse_cnt[end];
    pv->v = start;
    pv->w = weight;
    AK_reverse_cnt[end]++;

#ifdef TEST
    if (out_deg[start] > MAXD || in_deg[end] > MAXD) {
      cout << "exceed" << endl;
    }
#endif
  }

  //    cout << weights[9830][9481] << endl;

  //    n = unique(all_nodes, all_nodes + n) - all_nodes;

  TopoSortRemove();

  Vertex *g = AK, *r = AK_reverse;
  for (int i = 0; i < n; ++i) {
    sort(g, g + AK_cnt[i]);
    sort(r, r + AK_reverse_cnt[i]);
    g += MAXD;
    r += MAXD;
  }

#ifdef TEST
  auto end_time = chrono::steady_clock::now();
  cout << "读取并建立邻接表："
       << chrono::duration_cast<chrono::milliseconds>(end_time - start_time)
              .count()
       << " ms" << endl;
#endif
}

inline bool Check(int& x, int& y) {
  //    float ratio = y / (float)x;
  //    return ratio >= MIN_RATIO && ratio <= MAX_RATIO;
  //    return y <= 3 * x && x <= 5 * y;
  return y <= 3LL * x && x <= 5LL * y;
}

void WorkPart(int* res3, int* res4, int* res5, int* res6, int* res7,
              int res_cnt[4][5], int thread_no, Vertex* AK, uint8_t* AK_cnt,
              Vertex* AK_reverse, uint8_t* AK_reverse_cnt) {
#ifdef TEST
  auto start_time = chrono::steady_clock::now();
#endif
#ifdef TEST
  struct timeval tim {};
  gettimeofday(&tim, nullptr);
  double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
#endif

  int *p3 = res3, *p4 = res4, *p5 = res5, *p6 = res6, *p7 = res7;
  int v0, v1, v2, v3, v4, v5, v6;
  int w0_1, w1_2, w2_3, w3_4, w4_5, w5_6, w6_0;
  Vertex *g1, *g2, *g3, *g4, *g5, *g6, *g7;
  //    for (v0 = N - thread_no; v0 >= 0; v0 -= 4) {

  uint8_t* min_path = new uint8_t[n];
  int* last_ws = new int[n];

  for (v0 = n - 1; v0 >= 0; --v0) {
    if (v0 % 4 != thread_no) continue;
    if (AK_cnt[v0] == 0) continue;

    for (int i = v0 + 1; i < n; ++i) min_path[i] = 10;

    g6 = AK_reverse + v0 * MAXD;
    for (int i1 = 0; i1 < AK_reverse_cnt[v0]; ++i1) {
      v6 = g6->v;
      if (v6 <= v0) break;

      w6_0 = g6->w;
      last_ws[v6] = w6_0;

      min_path[v6] = 1;
      g5 = AK_reverse + v6 * MAXD;
      for (int i2 = 0; i2 < AK_reverse_cnt[v6]; ++i2) {
        v5 = g5->v;
        if (v5 <= v0) break;

        w5_6 = g5->w;
        if (Check(w5_6, w6_0)) {
          if (min_path[v5] > 2) min_path[v5] = 2;
          g4 = AK_reverse + v5 * MAXD;
          for (int i3 = 0; i3 < AK_reverse_cnt[v5]; ++i3) {
            v4 = g4->v;
            if (v4 <= v0) break;
            if (min_path[v4] > 3) min_path[v4] = 3;
            ++g4;
          }
        }
        ++g5;
      }

      ++g6;
    }

    // search for circuits
    g1 = AK + v0 * MAXD;
    for (int i1 = 0; i1 < AK_cnt[v0]; ++i1) {
      v1 = g1->v;
      if (v1 <= v0) break;

      w0_1 = g1->w;
      g2 = AK + v1 * MAXD;
      for (int i2 = 0; i2 < AK_cnt[v1]; ++i2) {
        v2 = g2->v;
        if (v2 <= v0) break;

        w1_2 = g2->w;
        if (v2 > v0 && Check(w0_1, w1_2)) {
          g3 = AK + v2 * MAXD;
          for (int i3 = 0; i3 < AK_cnt[v2]; ++i3) {
            v3 = g3->v;
            if (v3 < v0) break;

            w2_3 = g3->w;
            bool flag = false;
            if (v3 > v0 && v3 ^ v1 && (flag = Check(w1_2, w2_3))) {
              g4 = AK + v3 * MAXD;
              for (int i4 = 0; i4 < AK_cnt[v3]; ++i4) {
                v4 = g4->v;
                if (v4 < v0) break;

                w3_4 = g4->w;

                flag = false;

                if (v4 > v0 && min_path[v4] <= 3 && v4 ^ v2 && v4 ^ v1 &&
                    (flag = Check(w2_3, w3_4))) {
                  g5 = AK + v4 * MAXD;
                  for (int i5 = 0; i5 < AK_cnt[v4]; ++i5) {
                    v5 = g5->v;
                    if (v5 < v0) break;
                    w4_5 = g5->w;

                    flag = false;
                    if (v5 > v0 && min_path[v5] <= 2 && v5 ^ v3 && v5 ^ v2 &&
                        v5 ^ v1 && (flag = Check(w3_4, w4_5))) {
                      g6 = AK + v5 * MAXD;
                      for (int i6 = 0; i6 < AK_cnt[v5]; ++i6) {
                        v6 = g6->v;
                        if (v6 < v0) break;

                        w5_6 = g6->w;

                        flag = false;

                        if (v6 > v0 && min_path[v6] <= 1 && v6 ^ v4 &&
                            v6 ^ v3 && v6 ^ v2 && v6 ^ v1 &&
                            (flag = Check(w4_5, w5_6))) {
                          w6_0 = last_ws[v6];

                          if (Check(w5_6, w6_0) && Check(w6_0, w0_1)) {
                            p7[0] = v0;
                            p7[1] = v1;
                            p7[2] = v2;
                            p7[3] = v3;
                            p7[4] = v4;
                            p7[5] = v5;
                            p7[6] = v6;
                            p7 += 7;
                          }

                        } else if (v6 == v0 && Check(w5_6, w0_1) &&
                                   (flag || Check(w4_5, w5_6))) {
                          p6[0] = v0;
                          p6[1] = v1;
                          p6[2] = v2;
                          p6[3] = v3;
                          p6[4] = v4;
                          p6[5] = v5;
                          p6 += 6;
                          break;
                        }

                        ++g6;
                      }
                    } else if (v5 == v0 && Check(w4_5, w0_1) &&
                               (flag || Check(w3_4, w4_5))) {
                      p5[0] = v0;
                      p5[1] = v1;
                      p5[2] = v2;
                      p5[3] = v3;
                      p5[4] = v4;
                      p5 += 5;
                      break;
                    }

                    ++g5;
                  }
                } else if (v4 == v0 && Check(w3_4, w0_1) &&
                           (flag || Check(w2_3, w3_4))) {
                  p4[0] = v0;
                  p4[1] = v1;
                  p4[2] = v2;
                  p4[3] = v3;
                  p4 += 4;
                  break;
                }

                ++g4;
              }
            } else if (v3 == v0 && Check(w2_3, w0_1) &&
                       (flag || Check(w1_2, w2_3))) {
              p3[0] = v0;
              p3[1] = v1;
              p3[2] = v2;
              p3 += 3;
              break;
            }
            ++g3;
          }
        }
        ++g2;
      }

      ++g1;
    }
  }

  res_cnt[thread_no][0] = (p3 - res3) / 3;
  res_cnt[thread_no][1] = (p4 - res4) / 4;
  res_cnt[thread_no][2] = (p5 - res5) / 5;
  res_cnt[thread_no][3] = (p6 - res6) / 6;
  res_cnt[thread_no][4] = (p7 - res7) / 7;
#ifdef TEST
  auto end_time = chrono::steady_clock::now();
  cout << "Work: [ " << thread_no << " ]"
       << chrono::duration_cast<chrono::milliseconds>(end_time - start_time)
              .count()
       << " ms" << endl;
#endif
#ifdef LOCAL
  gettimeofday(&tim, nullptr);
  double t4 = tim.tv_sec + (tim.tv_usec / 1000000.0);
  printf("@ Work %d:\t[cost: %.4fs]\n", thread_no, t4 - t1);
#endif
}

void Work() {
#ifdef TEST
  auto start_time = chrono::steady_clock::now();
#endif

  thread t1(WorkPart, res3[0], res4[0], res5[0], res6[0], res7[0], res_cnt, 0,
            AK, AK_cnt, AK_reverse, AK_reverse_cnt);
  thread t2(WorkPart, res3[1], res4[1], res5[1], res6[1], res7[1], res_cnt, 1,
            AK, AK_cnt, AK_reverse, AK_reverse_cnt);
  thread t3(WorkPart, res3[2], res4[2], res5[2], res6[2], res7[2], res_cnt, 2,
            AK, AK_cnt, AK_reverse, AK_reverse_cnt);
  thread t4(WorkPart, res3[3], res4[3], res5[3], res6[3], res7[3], res_cnt, 3,
            AK, AK_cnt, AK_reverse, AK_reverse_cnt);

  t1.join();
  t2.join();
  t3.join();
  t4.join();

  //    WorkPart(res3[0], res4[0], res5[0], res6[0], res7[0], res_cnt, 0, AK,
  //    AK_cnt, AK_reverse, AK_reverse_cnt);

#ifdef TEST
  auto end_time = chrono::steady_clock::now();
  cout << "Work: "
       << chrono::duration_cast<chrono::milliseconds>(end_time - start_time)
              .count()
       << " ms" << endl;
#endif
}

void ResArrayToStr(int** res, char*& ans_cur, int len) {
  char* v;
  int* cur;
  int* cur1 = res[0] + (res_cnt[0][len - 3] - 1) * len;
  int* cur2 = res[1] + (res_cnt[1][len - 3] - 1) * len;
  int* cur3 = res[2] + (res_cnt[2][len - 3] - 1) * len;
  int* cur4 = res[3] + (res_cnt[3][len - 3] - 1) * len;

  int mn;

  while (cur1 >= res[0] || cur2 >= res[1] || cur3 >= res[2] || cur4 >= res[3]) {
    mn = INT_MAX;
    if (cur1 >= res[0] && mn > cur1[0]) mn = cur1[0];
    if (cur2 >= res[1] && mn > cur2[0]) mn = cur2[0];
    if (cur3 >= res[2] && mn > cur3[0]) mn = cur3[0];
    if (cur4 >= res[3] && mn > cur4[0]) mn = cur4[0];

    if (mn == cur1[0]) {
      cur = cur1;
      cur1 -= len;
    } else if (mn == cur2[0]) {
      cur = cur2;
      cur2 -= len;
    } else if (mn == cur3[0]) {
      cur = cur3;
      cur3 -= len;
    } else if (mn == cur4[0]) {
      cur = cur4;
      cur4 -= len;
    }

    for (int i = 0; i < len - 1; i++) {
      v = map_int_str[cur[i]];
      memcpy(ans_cur, v, map_int_len[cur[i]]);
      ans_cur += map_int_len[cur[i]];
    }
    v = map_int_str[cur[len - 1]];
    memcpy(ans_cur, v, map_int_len[cur[len - 1]]);
    ans_cur += map_int_len[cur[len - 1]];
    *(ans_cur - 1) = '\n';
  }
}

void SaveResultByFwrite(const char* outputFile) {
#ifdef TEST
  auto start_time = chrono::steady_clock::now();
#endif

  FILE* file = fopen(outputFile, "w");
  char* ans_cur = answer;
  char* ans_start = answer;

  int res_total_cnt = 0;
  for (int i = 0; i < 4; i++) {
    res_total_cnt += res_cnt[i][0] + res_cnt[i][1] + res_cnt[i][2] +
                     res_cnt[i][3] + res_cnt[i][4];
  }

#ifdef TEST
  cout << "总环数: " << res_total_cnt << endl;
#endif
  fprintf(file, "%d\n", res_total_cnt);
  //    fwrite("\n", 1, 1, file);

  ResArrayToStr(res3, ans_cur, 3);
  ResArrayToStr(res4, ans_cur, 4);
  ResArrayToStr(res5, ans_cur, 5);
  ResArrayToStr(res6, ans_cur, 6);
  ResArrayToStr(res7, ans_cur, 7);

  while (ans_start + 2097152 < ans_cur) {
    fwrite(ans_start, 1, 2097152, file);
    ans_start += 2097152;
  }

  fwrite(ans_start, 1, ans_cur - ans_start, file);
#ifdef TEST
  auto end_time = chrono::steady_clock::now();
  cout << "fwrite 写入结果："
       << chrono::duration_cast<chrono::milliseconds>(end_time - start_time)
              .count()
       << " ms" << endl;
#endif
}

int main() {
#ifdef TEST
  auto start_time = chrono::steady_clock::now();
  //    const char * test_file = "../data/test_data0.txt";
  //    const char * res_file = "../data/result0.txt";
  //    const char * test_file = "../data/test_data1.txt";
  //    const char * res_file = "../data/result1.txt";
  //    const char * test_file = "../data/test_data(1).txt";
  //    const char * res_file = "../data/result(1).txt";
  //    const char * test_file = "../data/test_data_2896262.txt";
  //    const char * res_file = "../data/result_2896262.txt";
  //    const char * test_file = "../data/test_data_300w.txt";
  //    const char * res_file = "../data/result_300w.txt";
  const char* test_file = "../data/19630345/test_data.txt";
  const char* res_file = "../data/19630345/result.txt";
//    const char * test_file = "../data/test_data_19882639.txt";
//    const char * res_file = "../data/result_19882639.txt";
#else
  const char* test_file = "/data/test_data.txt";
  const char* res_file = "/projects/student/result.txt";
#endif

  MyLoadData(test_file);

  Work();

  SaveResultByFwrite(res_file);

  //    sleep(5);

#ifdef TEST
  auto end_time = chrono::steady_clock::now();
  cout << "总耗时"
       << chrono::duration_cast<chrono::milliseconds>(end_time - start_time)
              .count()
       << " ms" << endl;
#endif
  exit(0);
}