#include <algorithm>
#include <cstdio>
#include <cstring>
#include <iostream>
#define mp make_pair
#include <queue>
using namespace std;
#define f2 double
#define f3 long double
#define rev(i, p) for (Edge *i = head[p]; i; i = i->next)
const int N = 1100;
struct Edge {
  int to, len;
  f3 val;
  Edge *next;
} * head[N], e[N << 3];
inline void addEdge(int x, int y, int z, f2 w) {
  static int _;
  Edge *t;
  t = &e[++_];
  t->to = y, t->len = z, t->val = w, t->next = head[x], head[x] = t;
  t = &e[++_];
  t->to = x, t->len = z, t->val = w, t->next = head[y], head[y] = t;
}
int a[N], ind[N], n, m;
int q[N], l, r, dis[N];
f2 ans[N];
f3 g[N], f[N];
priority_queue<pair<int, int> > pq;
bool vis[N];
void getTopp(int S) {
  q[r = 1] = S, l = 0;
  f[S] = 1;
  for (int i = 1; i <= n; i++) rev(j, i) {
      if (dis[i] + j->len == dis[j->to]) ind[j->to]++;
    }
  while (r > l) {
    int u = q[++l];
    rev(i, u) if (dis[u] + i->len == dis[i->to]) {
      if (!--ind[i->to]) q[++r] = i->to;
      f[i->to] += f[u] * i->val;
    }
  }
}
void getDis(int S) {
  dis[S] = 0;
  pq.push(mp(0, S));
  while (!pq.empty()) {
    int u = pq.top().second;
    pq.pop();
    if (vis[u]) continue;
    vis[u] = 1;
    rev(i, u) if (dis[i->to] > dis[u] + i->len) {
      dis[i->to] = dis[u] + i->len;
      pq.push(mp(-dis[i->to], i->to));
    }
  }
}
void getAns(int S) {
  for (int p = r; p > 1; p--) {
    int u = q[p];
    rev(i, u) if (dis[i->to] == dis[u] + i->len) {
      int v = i->to;
      g[u] += i->val * g[v];
    }
    ans[u] += a[S] * f[u] * g[u];
    g[u] += (f3)a[u] / f[u];
  }
}
int main() {
  scanf("%d%d", &n, &m);
  for (int i = 1; i <= n; i++) scanf("%d", &a[i]);
  for (int i = 1; i <= m; i++) {
    int x, y, z;
    double w;
    scanf("%d%d%d%lf", &x, &y, &z, &w);
    addEdge(x, y, z, w);
  }
  for (int i = 1; i <= n; i++) {
    memset(f, 0, sizeof f);
    memset(g, 0, sizeof g);
    memset(vis, 0, sizeof vis);
    memset(dis, 0x3f, sizeof dis);
    getDis(i), getTopp(i), getAns(i);
  }
  for (int i = 1; i <= n; i++) printf("%.8f\n", ans[i]);
}
#include <cstdio>#include <cstring>#include <iostream>#include <algorithm>#define mp make_pair#include <queue>using namespace std;
#define f2 double #define f3 long double #define rev(
    i, p) for (Edge *i = head[p]; i; i = i->next) const int N = 1100;
struct Edge {
  int to, len;
  f3 val;
  Edge *next;
} * head[N], e[N << 3];
inline void addEdge(int x, int y, int z, f2 w) {
  static int _;
  Edge *t;
  t = &e[++_];
  t->to = y, t->len = z, t->val = w, t->next = head[x], head[x] = t;
  t = &e[++_];
  t->to = x, t->len = z, t->val = w, t->next = head[y], head[y] = t;
}
int a[N], ind[N], n, m;
int q[N], l, r, dis[N];
f2 ans[N];
f3 g[N], f[N];
priority_queue<pair<int, int> > pq;
bool vis[N];
void getTopp(int S) {
  q[r = 1] = S, l = 0;
  f[S] = 1;
  for (int i = 1; i <= n; i++) rev(j, i) {
      if (dis[i] + j->len == dis[j->to]) ind[j->to]++;
    }
  while (r > l) {
    int u = q[++l];
    rev(i, u) if (dis[u] + i->len == dis[i->to]) {
      if (!--ind[i->to]) q[++r] = i->to;
      f[i->to] += f[u] * i->val;
    }
  }
}
void getDis(int S) {
  dis[S] = 0;
  pq.push(mp(0, S));
  while (!pq.empty()) {
    int u = pq.top().second;
    pq.pop();
    if (vis[u]) continue;
    vis[u] = 1;
    rev(i, u) if (dis[i->to] > dis[u] + i->len) {
      dis[i->to] = dis[u] + i->len;
      pq.push(mp(-dis[i->to], i->to));
    }
  }
}
void getAns(int S) {
  for (int p = r; p > 1; p--) {
    int u = q[p];
    rev(i, u) if (dis[i->to] == dis[u] + i->len) {
      int v = i->to;
      g[u] += i->val * g[v];
    }
    ans[u] += a[S] * f[u] * g[u];
    g[u] += (f3)a[u] / f[u];
  }
}
int main() {
  scanf("%d%d", &n, &m);
  for (int i = 1; i <= n; i++) scanf("%d", &a[i]);
  for (int i = 1; i <= m; i++) {
    int x, y, z;
    double w;
    scanf("%d%d%d%lf", &x, &y, &z, &w);
    addEdge(x, y, z, w);
  }
  for (int i = 1; i <= n; i++) {
    memset(f, 0, sizeof f);
    memset(g, 0, sizeof g);
    memset(vis, 0, sizeof vis);
    memset(dis, 0x3f, sizeof dis);
    getDis(i), getTopp(i), getAns(i);
  }
  for (int i = 1; i <= n; i++) printf("%.8f\n", ans[i]);
}