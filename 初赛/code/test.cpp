#include <string.h>

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
using namespace std;
struct node {
  int v, next;
} edge[280007];
int DFN[70000], LOW[70000];
int m_stack[70000], heads[70000], visit[70000], cnt, tot, m_index;
void add(int x, int y) {
  edge[++cnt].next = heads[x];
  edge[cnt].v = y;
  heads[x] = cnt;
  return;
}
int TOL = 0;
void tarjan(int x)  //代表第几个点在处理。递归的是点。
{
  DFN[x] = LOW[x] = ++tot;  // 新进点的初始化。
  m_stack[++m_index] = x;   //进站
  visit[x] = 1;             //表示在栈里
  for (int i = heads[x]; i != -1; i = edge[i].next) {
    if (!DFN[edge[i].v]) {  //如果没访问过
      tarjan(edge[i].v);    //往下进行延伸，开始递归
      LOW[x] = min(
          LOW[x],
          LOW[edge[i]
                  .v]);  //递归出来，比较谁是谁的儿子／父亲，就是树的对应关系，涉及到强连通分量子树最小根的事情。
    } else if (visit[edge[i].v]) {  //如果访问过，并且还在栈里。
      LOW[x] = min(LOW[x],
                   DFN[edge[i].v]);  //比较谁是谁的儿子／父亲。就是链接对应关系
    }
  }
  if (LOW[x] == DFN[x])  //发现是整个强连通分量子树里的最小根。
  {
    int ct = 0;
    do {
      // printf("%d ", m_stack[m_index]);
      visit[m_stack[m_index]] = 0;
      m_index--;
      ++ct;
    } while (x != m_stack[m_index + 1]);  //出栈，并且输出。
    // printf("\n");
    if (ct >= 1) ++TOL;
  }
  return;
}
int main() {
  memset(heads, -1, sizeof(heads));
  std::ifstream fin("../data/gen1/test_data.txt");
  int u, v, w;
  char c;
  while (fin >> u >> c >> v >> c >> w) {
    add(u, v);
  }
  fin.close();
  std::cerr << cnt << "\n";

  for (int i = 0; i < 70000; i++)
    if (!DFN[i] && heads[i] != -1)
      tarjan(i);  //当这个点没有访问过，就从此点开始。防止图没走完
  std::cerr << TOL << "\n";
  return 0;
}