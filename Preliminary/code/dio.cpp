#include <bits/stdc++.h>
#include <unistd.h>
using namespace std;

const string testCaseName = "1004812";

const int NODES = 280000 * 2;

// 静态图，第 i 个节点可以走到 pnt[fst[i]:lst[i]] 中的所有点
int fst[2][NODES], lst[2][NODES], pnt[2][NODES];

vector<vector<int>> answer[8];  // 答案

namespace Util {
double t_st = clock(), t_now = t_st;

void timing(string info = "no info") {
  cerr << "curr time: " << (clock() - t_st) / CLOCKS_PER_SEC << "s, "
       << "cost time: " << (clock() - t_now) / CLOCKS_PER_SEC << "s, " << info
       << "\n";
  t_now = clock();
}
}  // namespace Util

namespace IO {
void redirectToFile(string testid) {
  string inputPath = "/data/test_data.txt";
  string outputPath = "/projects/student/result.txt";

#ifdef _LOCAL_
  cerr << "testing case: " + testid << "\n";
  inputPath = "../data/" + testid + "/test_data.txt";
  outputPath = "../data/" + testid + "/myresult.txt";
#endif

  freopen(inputPath.c_str(), "r", stdin);
  freopen(outputPath.c_str(), "w", stdout);
}
pair<int, int>
    rd[NODES / 2];  // rd = record，转账记录：本端账号，对端账号，金额（暂无）

int posInReadBuf[NODES], numberLength[NODES];

char readBuf[NODES * 18];
int rptr = 0;
pair<int, int> readIn() {
  readBuf[fread(readBuf, 1, NODES * 18, stdin)] = 0;

  int rds = 0, dcs = 0;

  for (int x; readBuf[rptr]; rds++) {
    x = 0;
    while (readBuf[rptr] != ',') {
      x = x * 10 + readBuf[rptr++] - '0';
    }
    rd[rds].first = x;
    posInReadBuf[x] = rptr;

    ++rptr;

    x = 0;
    while (readBuf[rptr] != ',') {
      x = x * 10 + readBuf[rptr++] - '0';
    }
    rd[rds].second = x;
    posInReadBuf[x] = rptr;

    while (readBuf[rptr++] != '\n')
      ;

    dcs = max(dcs, max(rd[rds].first, rd[rds].second));
  }

  for (int i = 1; i <= dcs; ++i) {
    numberLength[i] = numberLength[i / 10] + 1;
  }
  numberLength[0] = 1;
  for (int i = 0; i <= dcs; ++i) {
    posInReadBuf[i] -= numberLength[i];
  }

  return {rds, dcs + 1};
}

void buildGraph(int n, int m) {
  // TODO: 这里的两个排序改成基数排序
  sort(rd, rd + m, [](const pair<int, int> &a, const pair<int, int> &b) {
    return a.first == b.first ? a.second < b.second : a.first < b.first;
  });
  for (int i = 0, eid = 0; i < n; ++i) {
    fst[0][i] = eid;
    while (eid < m && rd[eid].first == i) {
      pnt[0][eid] = rd[eid].second;
      ++eid;
    }
    lst[0][i] = eid;
  }

  sort(rd, rd + m, [](const pair<int, int> &a, const pair<int, int> &b) {
    return a.second == b.second ? a.first < b.first : a.second < b.second;
  });
  for (int i = 0, eid = 0; i < n; ++i) {
    fst[1][i] = eid;
    while (eid < m && rd[eid].second == i) {
      pnt[1][eid] = rd[eid].first;
      ++eid;
    }
    lst[1][i] = eid;
  }
}

char writeBuf[3000000 * 80];
int wptr = 0;  // wptr 表示下一个可以插入的位置
inline void writeInt(int t) {
  memcpy(writeBuf + wptr, readBuf + posInReadBuf[t], numberLength[t]);
  wptr += numberLength[t];
}
void printAnswer() {
  int tot = 0;
  for (int length = 3; length <= 7; ++length) {
    tot += answer[length].size();
  }

  cerr << "total loops = " << tot << "\n";

  wptr = sprintf(writeBuf, "%d\n", tot);
  for (int length = 3; length <= 7; ++length) {
    // vi 代表每个环
    for (const auto &vi : answer[length]) {
      for (int x : vi) {
        writeInt(x);
        writeBuf[wptr++] = ',';
      }
      writeBuf[wptr - 1] = '\n';
    }
  }

  Util::timing("wrtieBuf prepare complete.");

  fwrite(writeBuf, 1, wptr, stdout);
}

}  // namespace IO

namespace Solve {
int curr;
int disFrom[4][NODES];  // disFrom[i][j]=curr 表示到 j 到 curr 存在距离等于 i+1
                        // 的路径

int que[4][NODES], qp[4];

void bfs() {
  for (int i = 0; i < 4; ++i) {
    qp[i] = 0;
  }
  que[0][qp[0]++] = curr;
  for (int dis = 0; dis < 4; ++dis) {
    while (qp[dis]) {
      int u = que[dis][--qp[dis]];
      for (int &eid = fst[1][u]; eid < lst[1][u] && pnt[1][eid] <= curr; ++eid)
        ;
      for (int eid = fst[1][u]; eid < lst[1][u]; ++eid) {
        if (disFrom[dis][pnt[1][eid]] != curr) {
          disFrom[dis][pnt[1][eid]] = curr;
          if (dis < 3) {
            que[dis + 1][qp[dis + 1]++] = pnt[1][eid];
          }
        }
      }
    }
  }
}

bool occ[NODES];
vector<int> tmp;
void dfs(int u, int step) {
  /*
  目标：找到长度 3-7 的环
  结果判断条件：step 为2-6时，看一下有没有通向起点的点 (disFrom 1 是否为 true)
  边界：step为6时返回
  剪枝：
  1. step 为6时，如果 disFrom1 存在就记录然后返回，否则就返回。
  2. step 为5时，如果 disFrom1 存在就记录，如果disFrom2存在就递归，否则返回
  3. step 为4时，如果 disFrom1 存在就记录，如果disFrom23存在就递归，否则返回
  4. step 为3时，如果 disFrom1 存在就记录，如果disfrom234存在就递归，否则返回
  5. step 为2时，如果 disfrom1 存在就记录，然后必须递归
  6. step 为0或1时，直接递归
  */

  tmp.push_back(u);
  occ[u] = 1;

  bool needRecur = 1;

  // 记录环

  if (step >= 2) {
    if (disFrom[0][u] == curr) {
      answer[step + 1].push_back(tmp);
    }
    if (step == 6) {
      needRecur = 0;
    } else if (step == 5) {
      needRecur = disFrom[1][u] == curr;
    } else if (step == 4) {
      needRecur = disFrom[2][u] == curr || disFrom[1][u] == curr;
    } else if (step == 3) {
      needRecur = disFrom[3][u] == curr || disFrom[2][u] == curr ||
                  disFrom[1][u] == curr;
    }
  }

  if (needRecur) {
    for (int &eid = fst[0][u]; eid < lst[0][u] && pnt[0][eid] <= curr; ++eid)
      ;
    for (int eid = fst[0][u]; eid < lst[0][u]; ++eid) {
      if (occ[pnt[0][eid]] == 0) {
        dfs(pnt[0][eid], step + 1);
      }
    }
  }

  occ[u] = 0;
  tmp.pop_back();
}

// 找到包含点 s，且节点编号均大于等于 s 的所有环，按长度存入 answer 中
void solve(int s) {
  // 剪枝入度或出度为0的点
  if (lst[0][s] == fst[0][s] || lst[1][s] == fst[1][s]) {
    return;
  }

  curr = s;

  bfs();
  dfs(s, 0);
}
}  // namespace Solve

int main(void) {
  IO::redirectToFile(testCaseName);  // 重定向到文件

  int n, m;
  tie(m, n) = IO::readIn();  // 读入并离散化，返回边数（记录数），点数（账号数）

  cerr << "n, m = " << n << ", " << m << "\n";

  IO::buildGraph(n, m);  // 构建图

  Util::timing("build graph suc.");

  memset(Solve::disFrom, -1, sizeof(Solve::disFrom));
  for (int i = 0; i < n; ++i) {
    Solve::solve(i);
  }

  Util::timing("find loop suc.");

  // 对结果排序并去重
  // for(int length=3; length<=7; ++length) {
  // 	auto &v = answer[length];
  // 	sort(v.begin(), v.end());
  // 	v.resize(unique(v.begin(),v.end())-v.begin());
  // }

  //输出
  IO::printAnswer();

  Util::timing("output suc.");

  return 0;
}