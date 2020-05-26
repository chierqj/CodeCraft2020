#include <algorithm>
#include <cstdio>
#include <cstring>

using namespace std;

#define ll long long
#define MAXN 1000001

struct edge {
  edge *next;
  int t, d;
} * head[MAXN];

void AddEdge(int s, int t, int d) {
  edge *p = new (edge);
  p->t = t, p->d = d, p->next = head[s];
  head[s] = p;
}

ll dist[MAXN];
bool f[MAXN];
const ll inf = (ll)(0x7fffffff) * (ll)(0x7fffffff);
int n, m;
ll T, rxa, rxc, rya, ryc, rp;

struct node {
  int left, right, child;
  node() { left = right = child; }
} h[MAXN];

int roof = 0;

int Join(int v, int u) {
  if (dist[v] < dist[u]) swap(v, u);
  h[v].left = u, h[v].right = h[u].child, h[h[u].child].left = v;
  h[u].child = v;
  return u;
}

void Push(int v) {
  if (!roof)
    roof = v;
  else
    roof = Join(roof, v);
}

int Top() { return roof; }

void Update(int v) {
  if (v != roof) {
    if (h[h[v].left].child == v) {
      h[h[v].left].child = h[v].right;
    } else {
      h[h[v].left].right = h[v].right;
    }
    if (h[v].right) h[h[v].right].left = h[v].left;
    h[v].left = h[v].right = 0;
    roof = Join(roof, v);
  }
}

int sta[MAXN], top;

void Pop() {
  if (!h[roof].child)
    roof = 0;
  else {
    top = 0;
    int t = h[roof].child;
    while (t) {
      if (h[t].right) {
        int k = h[h[t].right].right;
        int v = h[t].right;
        h[t].left = h[t].right = h[v].left = h[v].right = 0;
        sta[++top] = Join(v, t);
        t = k;
      } else {
        sta[++top] = t;
        h[t].left = h[t].right = 0;
        break;
      }
    }
    roof = sta[top];
    for (int i = top - 1; i; --i) roof = Join(roof, sta[i]);
  }
}

void Dijstra() {
  memset(f, false, sizeof(f));
  for (int i = 0; i++ < n;) dist[i] = inf;
  dist[1] = 0, Push(1), f[1] = true;
  for (int i = 0; i++ < n - 1;) {
    int v = Top();
    Pop(), f[v] = false;
    if (v == n) break;
    for (edge *p = head[v]; p; p = p->next)
      if (dist[p->t] > dist[v] + (ll)(p->d)) {
        dist[p->t] = dist[v] + (ll)(p->d);
        if (!f[p->t])
          Push(p->t), f[p->t] = true;
        else
          Update(p->t);
      }
  }
}

int main() {
  scanf("%d%d", &n, &m);
  scanf("%lld%lld%lld%lld%lld%lld", &T, &rxa, &rxc, &rya, &ryc, &rp);
  memset(head, 0, sizeof(head));
  ll x = 0, y = 0, z = 0;
  for (int i = 0; i++ < T;) {
    x = (x * rxa + rxc) % rp;
    y = (y * rya + ryc) % rp;
    ll a = x % n + 1, b = y % n + 1;
    if (a > b) swap(a, b);
    ll d = (ll)(100000000) - 100 * a;
    AddEdge(a, b, d);
  }
  for (int i = 0; i++ < m - T;) {
    int s, t;
    ll d;
    scanf("%d%d%lld", &s, &t, &d);
    AddEdge(s, t, d);
  }
  Dijstra();
  printf("%lld\n", dist[n]);
  return 0;
}