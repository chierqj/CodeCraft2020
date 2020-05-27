#define FASTER
#ifdef FASTER
#pragma GCC diagnostic error "-std=c++11"
#pragma GCC optimize(2)
#pragma GCC optimize(3)
#pragma GCC optimize("Ofast")
#pragma GCC optimize("inline")
#pragma GCC optimize("-fgcse")
#pragma GCC optimize("-fgcse-lm")
#pragma GCC optimize("-fipa-sra")
#pragma GCC optimize("-ftree-pre")
#pragma GCC optimize("-ftree-vrp")
#pragma GCC optimize("-fpeephole2")
#pragma GCC optimize("-ffast-math")
#pragma GCC optimize("-fsched-spec")
#pragma GCC optimize("unroll-loops")
#pragma GCC optimize("-falign-jumps")
#pragma GCC optimize("-falign-loops")
#pragma GCC optimize("-falign-labels")
#pragma GCC optimize("-fdevirtualize")
#pragma GCC optimize("-fcaller-saves")
#pragma GCC optimize("-fcrossjumping")
#pragma GCC optimize("-fthread-jumps")
#pragma GCC optimize("-funroll-loops")
#pragma GCC optimize("-fwhole-program")
#pragma GCC optimize("-freorder-blocks")
#pragma GCC optimize("-fschedule-insns")
#pragma GCC optimize("inline-functions")
#pragma GCC optimize("-ftree-tail-merge")
#pragma GCC optimize("-fschedule-insns2")
#pragma GCC optimize("-fstrict-aliasing")
#pragma GCC optimize("-fstrict-overflow")
#pragma GCC optimize("-falign-functions")
#pragma GCC optimize("-fcse-skip-blocks")
#pragma GCC optimize("-fcse-follow-jumps")
#pragma GCC optimize("-fsched-interblock")
#pragma GCC optimize("-fpartial-inlining")
#pragma GCC optimize("no-stack-protector")
#pragma GCC optimize("-freorder-functions")
#pragma GCC optimize("-findirect-inlining")
#pragma GCC optimize("-fhoist-adjacent-loads")
#pragma GCC optimize("-frerun-cse-after-loop")
#pragma GCC optimize("inline-small-functions")
#pragma GCC optimize("-finline-small-functions")
#pragma GCC optimize("-ftree-switch-conversion")
#pragma GCC optimize("-foptimize-sibling-calls")
#pragma GCC optimize("-fexpensive-optimizations")
#pragma GCC optimize("-funsafe-loop-optimizations")
#pragma GCC optimize("inline-functions-called-once")
#pragma GCC optimize("-fdelete-null-pointer-checks")
#endif
#include <bits/extc++.h>
#include <bits/stdc++.h>
// luogu-judger-enable-o2
#include <bits/stdc++.h>

#define _rep(_1, _2, _3, _4, name, ...) name
#define rep2(i, n) rep3(i, 0, n)
#define rep3(i, a, b) rep4(i, a, b, 1)
#define rep4(i, a, b, c) for (int i = int(a); i < int(b); i += int(c))
#define rep(...) _rep(__VA_ARGS__, rep4, rep3, rep2, _)(__VA_ARGS__)

using namespace std;

using i64 = long long;
using u8 = unsigned char;
using u32 = unsigned;
using u64 = unsigned long long;
using f80 = long double;

#ifndef __APPLE__
#define fread fread_unlocked
#define fwrite fwrite_unlocked
#endif

class IO {
  enum { E = 3, N = 1000, SIZE = 1 << 14 };

 public:
  IO() : ii(SIZE), oi(0) { init(); }
  ~IO() {
    if (oi) flush(oi);
  }
  void init() {
    rep(i, N) {
      int n = i, l = 1;
      rep(j, E) num[i][E - 1 - j] = n % 10 + '0', n /= 10, l += n > 0;
      offsets[i] = E - l;
    }
  }
  char rchar() {
    if (ii == SIZE) {
      ii = 0;
      int s = fread(in, 1, SIZE, stdin);
      if (s < SIZE) fill(in + s, in + SIZE, EOF);
    }
    return in[ii++];
  }
  int read_int() {
    char c;
    int ret = 0, sign = 0;
    while ((c = rchar()) < '-') {
      if (c == EOF) return -1;
    }
    if (c == '-')
      sign = 1;
    else
      ret = c - '0';
    while ((c = rchar()) >= '0') ret = ret * 10 + c - '0';
    return sign ? -ret : ret;
  }
  i64 read_i64() {
    char c;
    i64 ret = 0;
    while ((c = rchar()) < '0')
      ;
    ret = c - '0';
    while ((c = rchar()) >= '0') ret = ret * 10 + c - '0';
    return ret;
  }
  u32 read_u32() {
    char c;
    u32 ret = 0;
    while ((c = rchar()) < '0')
      ;
    ret = c - '0';
    while ((c = rchar()) >= '0') ret = ret * 10 + c - '0';
    return ret;
  }
  void write_u32(u32 n, bool nl = true) {
    const u32 ten9 = 1e9;
    if (n >= ten9) {
      wchar('0' + n / ten9);
      n %= ten9;
      u32 q = n / 1000000;
      n %= 1000000;
      rep(i, 3) wchar(num[q][i]);
      q = n / 1000;
      n %= 1000;
      rep(i, 3) wchar(num[q][i]);
      rep(i, 3) wchar(num[n][i]);
    } else if (n >= u32(1e6)) {
      u32 q = n / 1000000;
      n %= 1000000;
      rep(i, offsets[q], 3) wchar(num[q][i]);
      q = n / 1000;
      n %= 1000;
      rep(i, 3) wchar(num[q][i]);
      rep(i, 3) wchar(num[n][i]);
    } else if (n >= u32(1e3)) {
      u32 q = n / 1000;
      n %= 1000;
      rep(i, offsets[q], 3) wchar(num[q][i]);
      rep(i, 3) wchar(num[n][i]);
    } else {
      rep(i, offsets[n], 3) wchar(num[n][i]);
    }
    if (nl) wchar('\n');
  }
  void write_int(int n, bool nl = true) {
    if (n < 0) {
      n = -n;
      wchar('-');
    }
    write_u32(n, nl);
  }
  void write_u32z(u32 n) {
    u32 q = n / 1000000;
    n %= 1000000;
    rep(i, 3) wchar(num[q][i]);
    q = n / 1000;
    n %= 1000;
    rep(i, 3) wchar(num[q][i]);
    rep(i, 3) wchar(num[n][i]);
  }
  void write_u64(u64 n, bool nl = true) {
    static const u64 t18 = u64(1e18);
    static const u32 t9 = u32(1e9);
    if (n >= t18) {
      wchar('0' + n / t18);
      n %= t18;
      write_u32z(n / t9);
      write_u32z(n % t9);
    } else if (n >= t9) {
      write_u32(n / t9, false);
      write_u32z(n % t9);
    } else {
      write_u32(n, false);
    }
    if (nl) wchar('\n');
  }
  void write_i64(i64 n, bool nl = true) {
    if (n < 0) {
      n = -n;
      wchar('-');
    }
    write_u64(n, nl);
  }
  void wchar(char c) {
    out[oi++] = c;
    if (oi == SIZE) flush(oi);
  }
  void wstr(const char* str) {
    while (*str) wchar(*str++);
  }
  void wstr(const char* str, int len) { rep(i, len) wchar(str[i]); }
  void flush(int size) {
    fwrite(out, 1, size, stdout);
    oi = 0;
  }
  int ii, oi;
  int offsets[N];
  char num[N][E], in[SIZE], out[SIZE];
} io;

template <typename I>
void umin(I& a, I b) {
  if (a > b) a = b;
}

const int inf = numeric_limits<int>::max();

struct RadixHeap {
  RadixHeap() : sz(0), last(0) { fill(vm, vm + 33, inf); }
  static constexpr int pos(int x) { return x == 0 ? 0 : 32 - __builtin_clz(x); }
  void push(int x, int y) {
    ++sz;
    const auto vi = pos(x ^ last);
    v[vi].emplace_back(x, y);
    umin(vm[vi], x);
  }
  pair<int, int> pop() {
    const auto ret = get();
    v[0].pop_back();
    --sz;
    return ret;
  }
  pair<int, int> top() { return get(); }
  pair<int, int> get() {
    if (v[0].empty()) {
      int i = 1;
      for (; v[i].empty(); ++i)
        ;
      last = vm[i];
      for (const auto& p : v[i]) {
        const auto vi = pos(p.first ^ last);
        v[vi].emplace_back(p);
        umin(vm[vi], p.first);
      }
      v[i].clear();
      vm[i] = inf;
    }
    return v[0].back();
  }
  bool empty() const { return sz == 0; }
  int size() const { return sz; }
  vector<pair<int, int> > v[33];
  int vm[33];
  int sz, last;
};

const int N_MAX = 100010;
const int M_MAX = 200010;

struct Edge {
  int from, to, cost;
};
int ofs[N_MAX];
pair<int, int> edges[M_MAX];
int dist[N_MAX];

void solve() {
  int N = io.read_u32(), M = io.read_u32(), s = io.read_u32() - 1;
  vector<Edge> in(M);
  rep(i, M) {
    int u = io.read_u32(), v = io.read_u32(), c = io.read_u32();
    ++ofs[u];
    in[i] = {u - 1, v - 1, c};
  }

  for (int i = 0; i < N; ++i) ofs[i + 1] += ofs[i];
  for (const auto& e : in) edges[ofs[e.from]++] = {e.to, e.cost};
  for (int i = N - 1; i > 0; --i) ofs[i] = ofs[i - 1];
  ofs[0] = 0;

  RadixHeap heap;
  heap.push(0, s);
  fill(dist, dist + N, inf);
  dist[s] = 0;
  while (!heap.empty()) {
    int v, d;
    tie(d, v) = heap.pop();
    if (d > dist[v]) continue;
    for (int ei = ofs[v]; ei < ofs[v + 1]; ++ei) {
      int u = edges[ei].first;
      int w = d + edges[ei].second;
      if (w < dist[u]) {
        dist[u] = w;
        heap.push(w, u);
      }
    }
  }
  rep(i, N) {
    io.write_u32(dist[i], false);
    io.wchar(i == N - 1 ? '\n' : ' ');
  }
}

int main() {
  solve();
  return 0;
}