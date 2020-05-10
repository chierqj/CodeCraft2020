#define _CRT_SECURE_NO_DEPRECATE
#include <time.h>

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <queue>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#define PIECE 64

using namespace std;

typedef unsigned int ui;
typedef long long unsigned int lli;

int path7[4][20000000][7] = {0};
int path6[4][2500000][6] = {0};
int path5[4][2500000][5] = {0};
int path4[4][500000][4] = {0};
int path3[4][500000][3] = {0};
int visited[4][500000] = {0};
char out[2560000];

struct Path {
  int l0;
  int l1;
  int l2;
  ui m;
  ui mo;
  Path() : l0(-1), l1(-1), l2(-1), m(0), mo(0) {}
};

struct Record {
  ui v;
  ui m;
  Record() : v(), m() {}
  bool operator<(const Record& rhs) const { return v < rhs.v; }
};

vector<vector<Path>> preCord[PIECE];

class Solution {
 public:
  vector<vector<Record>> Gfront;  //��һ�ν�ͼ�Ľ����������Ըĳ�����
  vector<vector<Record>> Graph;  //���˽���֮��ͼ��������Ըĳ�����
  vector<ui> ids;     //��IDתΪ�˻�
  vector<ui> inputs;  //�����ת�˼�¼
  vector<ui> inputs_m;
  vector<ui> inDegree;   //ÿ���ڵ�����
  vector<ui> outDegree;  //ÿ���ڵ�����
  vector<string> idsComma;  // 0...n to sorted id�����������ļ����
  vector<string> idsLF;  // 0...n to sorted id�����������ļ����
  vector<vector<vector<pair<int, int>>>> p4l;  //�洢����ת��3�ε�·��

  int path7count[4] = {0};  //��������
  int path6count[4] = {0};
  int path5count[4] = {0};
  int path4count[4] = {0};
  int path3count[4] = {0};
  int si;

  ui max_ID;

  int nodeCount;
  int tnCount;
  int pathCount[4] = {0};
  int Gcount;
  int reGcount;

  void parseInput(string& testFile) {
    FILE* file = fopen(testFile.c_str(), "rb");
    fseek(file, 0, SEEK_END);
    size_t size = ftell(file);

    fseek(file, 0, SEEK_SET);
    char* rec = new char[size];
    size_t cc = fread(rec, size, sizeof(char), file);
    ui u, v, c;
    inputs.reserve(4000000);
    inputs_m.reserve(2000000);
    int ss = 0;
    while (1) {
      u = atoi(rec + ss);
      while (rec[ss] != ',') {
        ss++;
        if (ss >= size - 1) {
          break;
        }
      }
      ss++;

      v = atoi(rec + ss);
      while (rec[ss] != ',') {
        ss++;
        if (ss >= size - 1) {
          break;
        }
      }
      ss++;

      c = atoi(rec + ss);
      while (rec[ss] != '\n') {
        ss++;
        if (ss >= size - 1) {
          break;
        }
      }
      ss++;

      inputs.push_back(u);
      inputs.push_back(v);
      inputs_m.push_back(c);

      if (ss >= size - 1) {
        break;
      }
    }
  }

  void constructGraph() {
    ids = inputs;
    int record_size = inputs.size();
    sort(ids.begin(), ids.end());
    ids.erase(unique(ids.begin(), ids.end()), ids.end());
    max_ID = *(ids.end() - 1);
    cout << max_ID << endl;
    cout << ids.size() << endl;
    if (max_ID < 400000000) {
      ui* idHash = new ui[max_ID + 1];
      nodeCount = 0;
      for (ui& x : ids) {
        idHash[x] = nodeCount++;
      }
      Gcount = 0;
      Gfront = vector<vector<Record>>(nodeCount);
      inDegree = vector<ui>(nodeCount, 0);
      outDegree = vector<ui>(nodeCount, 0);
      Record Rec;
      for (int i = 0; i < record_size; i += 2) {
        ui u = idHash[inputs[i]];
        ui v = idHash[inputs[i + 1]];
        ui m = inputs_m[i / 2];
        Rec.v = v;
        Rec.m = m;
        Gfront[u].push_back(Rec);
        ++inDegree[v];
        ++outDegree[u];
      }
      delete[] idHash;
    } else {
      unordered_map<ui, int> idHash;
      nodeCount = 0;
      for (ui& x : ids) {
        idHash[x] = nodeCount++;
      }
      Gcount = 0;
      Gfront = vector<vector<Record>>(nodeCount);
      inDegree = vector<ui>(nodeCount, 0);
      outDegree = vector<ui>(nodeCount, 0);
      Record Rec;
      for (int i = 0; i < record_size; i += 2) {
        ui u = idHash[inputs[i]];
        ui v = idHash[inputs[i + 1]];
        ui m = inputs_m[i / 2];
        Rec.v = v;
        Rec.m = m;
        Gfront[u].push_back(Rec);
        ++inDegree[v];
        ++outDegree[u];
      }
    }
  }

  void topologicalSort() {
    queue<int> topoQueue;
    int topocount = 0;
    for (int i = 0; i < nodeCount; i++) {
      if (inDegree[i] == 0) {
        topoQueue.push(i);
      }
    }
    while (!topoQueue.empty()) {
      int node = topoQueue.front();
      topoQueue.pop();
      for (auto& t : Gfront[node]) {
        inDegree[t.v]--;
        if (!inDegree[t.v]) {
          topoQueue.push(t.v);
        }
      }
    }
    for (int i = 0; i < nodeCount; i++) {
      if (inDegree[i] == 0) {
        Gfront[i].clear();
        topocount++;
      }
    }

    for (int i = 0; i < nodeCount; i++) {
      if (outDegree[i] == 0) {
        topoQueue.push(i);
      }
    }
    while (!topoQueue.empty()) {
      int node = topoQueue.front();
      topoQueue.pop();
      for (auto& t : Gfront[node]) {
        outDegree[t.v]--;
        if (!outDegree[t.v]) {
          topoQueue.push(t.v);
        }
      }
    }
    for (int i = 0; i < nodeCount; i++) {
      if (outDegree[i] == 0) {
        Gfront[i].clear();
        topocount++;
      }
    }
    cout << "topo delete " << topocount << endl;
  }

  void reConstructGraph() {
    if (max_ID < 400000000) {
      int* tidHash = new int[max_ID + 1];
      tnCount = 0;
      for (int i = 0; i < nodeCount; i++) {
        if (!Gfront[i].empty()) {
          tidHash[ids[i]] = tnCount++;
          idsComma.push_back(to_string(ids[i]) + ',');
          idsLF.push_back(to_string(ids[i]) + '\n');
        }
      }

      Graph = vector<vector<Record>>(tnCount);

      Record rec;

      for (int i = 0; i < nodeCount; i++) {
        if (!Gfront[i].empty()) {
          int u = tidHash[ids[i]];
          for (auto& t : Gfront[i]) {
            if (!Gfront[t.v].empty()) {
              rec.v = tidHash[ids[t.v]];
              rec.m = t.m;
              Graph[u].push_back(rec);
            }
          }
        }
      }
    } else {
      unordered_map<ui, int> tidHash;
      tnCount = 0;
      for (int i = 0; i < nodeCount; i++) {
        if (!Gfront[i].empty()) {
          tidHash[ids[i]] = tnCount++;
          idsComma.push_back(to_string(ids[i]) + ',');
          idsLF.push_back(to_string(ids[i]) + '\n');
        }
      }

      Graph = vector<vector<Record>>(tnCount);

      Record rec;

      for (int i = 0; i < nodeCount; i++) {
        if (!Gfront[i].empty()) {
          int u = tidHash[ids[i]];
          for (auto& t : Gfront[i]) {
            if (!Gfront[t.v].empty()) {
              rec.v = tidHash[ids[t.v]];
              rec.m = t.m;
              Graph[u].push_back(rec);
            }
          }
        }
      }
    }

    for (int i = 0; i < tnCount; i++) {
      sort(Graph[i].begin(), Graph[i].end());
    }
  }

  inline bool com(ui& X, ui& Y) {
    return ((lli)Y * 5 >= (lli)X) && ((lli)Y <= (lli)X * 3);
  }

  void ConstructPrecordTh(int ID) {
    for (int i = ID; i < tnCount; i += 4) {
      int t = i % PIECE;
      Path pt;
      pt.l0 = i;
      for (auto& level2 : Graph[i]) {
        pt.m = level2.m;
        if (i > level2.v) {
          pt.mo = level2.m;
          preCord[t][level2.v].push_back(pt);
        }
        for (auto& level3 : Graph[level2.v]) {
          if (com(level2.m, level3.m)) {
            if (i > level3.v && level2.v > level3.v) {
              {
                pt.l1 = level2.v;
                pt.mo = level3.m;
                preCord[t][level3.v].push_back(pt);
              }
            }
            for (auto& level4 : Graph[level3.v]) {
              if (i > level4.v && level2.v > level4.v && level3.v > level4.v) {
                if (com(level3.m, level4.m)) {
                  pt.l1 = level2.v;
                  pt.l2 = level3.v;
                  pt.mo = level4.m;
                  preCord[t][level4.v].push_back(pt);
                }
              } else
                break;
            }
            pt.l2 = -1;
          }
        }
        pt.l1 = -1;
      }
    }
  }

  void ConstructPrecord() {
    for (int i = 0; i < PIECE; i++) {
      preCord[i].resize(tnCount);
    }
    thread t0(&Solution::ConstructPrecordTh, this, 0);
    thread t1(&Solution::ConstructPrecordTh, this, 1);
    thread t2(&Solution::ConstructPrecordTh, this, 2);
    thread t3(&Solution::ConstructPrecordTh, this, 3);
    t0.join();
    t1.join();
    t2.join();
    t3.join();
  }

  void Find_Path(int* path, Record root, int level1, int id, ui m) {
    int t = root.v % PIECE;
    for (auto& x : preCord[t][level1]) {
      if (x.l0 == root.v && com(root.m, x.m) && com(x.mo, m)) {
        path[4] = root.v;
        if (x.l1 == -1) {
          pathCount[id]++;
          memcpy(path5[id][path5count[id]], path, 5 * si);
          path5count[id]++;
        } else if (!visited[id][x.l1]) {
          path[5] = x.l1;
          if (x.l2 == -1) {
            pathCount[id]++;
            memcpy(path6[id][path6count[id]], path, 6 * si);
            path6count[id]++;
          } else if (!visited[id][x.l2]) {
            path[6] = x.l2;
            pathCount[id]++;
            memcpy(path7[id][path7count[id]], path, 7 * si);
            path7count[id]++;
          }
        }
      }
      if (x.l0 > root.v) break;
    }
  }

  void DFS(int id) {
    int depth = 0;
    int path[7];
    for (int level1 = id; level1 < tnCount; level1 = level1 + 4) {
      // if (level1 % 500 == 0) {
      //   cout << level1 << " node has comleted" << endl;
      // }
      path[0] = level1;
      for (auto& level2 : Graph[level1]) {
        depth = 1;
        path[1] = level2.v;
        if (level2.v > level1) {
          visited[id][level2.v] = 1;
          for (auto& level3 : Graph[level2.v]) {
            if (level3.v > level1 && com(level2.m, level3.m)) {
              depth = 2;
              path[2] = level3.v;
              visited[id][level3.v] = 1;
              for (auto& level4 : Graph[level3.v]) {
                if (!visited[id][level4.v] && com(level3.m, level4.m)) {
                  depth = 3;
                  path[3] = level4.v;
                  visited[id][level4.v] = 1;
                  if (level4.v == level1 && com(level4.m, level2.m)) {
                    pathCount[id]++;
                    memcpy(path3[id][path3count[id]], path, 3 * si);
                    path3count[id]++;
                  }
                  if (level4.v > level1) {
                    for (auto& level5 : Graph[level4.v]) {
                      if (!visited[id][level5.v] && com(level4.m, level5.m)) {
                        depth = 4;
                        path[4] = level5.v;
                        visited[id][level5.v] = 1;
                        if (level5.v == level1 && com(level5.m, level2.m)) {
                          pathCount[id]++;
                          memcpy(path4[id][path4count[id]], path, 4 * si);
                          path4count[id]++;
                        }
                        if (level5.v > level1) {
                          Find_Path(path, level5, level1, id, level2.m);
                        }
                        visited[id][level5.v] = 0;
                      }
                    }
                  }
                  visited[id][level4.v] = 0;
                }
              }
              visited[id][level3.v] = 0;
            }
          }
          visited[id][level2.v] = 0;
        }
      }
    }
  }

  void solve() {
    si = sizeof(int);
    thread t0(&Solution::DFS, this, 0);
    thread t1(&Solution::DFS, this, 1);
    thread t2(&Solution::DFS, this, 2);
    thread t3(&Solution::DFS, this, 3);
    t0.join();
    t1.join();
    t2.join();
    t3.join();
    cout << pathCount[0] + pathCount[1] + pathCount[2] + pathCount[3] << endl;
  }

  void save_fwrite(const string& outputFile) {
    FILE* fp = fopen(outputFile.c_str(), "wb");
    char buf[1024];
    int pathTotal = pathCount[0] + pathCount[1] + pathCount[2] + pathCount[3];
    cout << "pathtotal " << pathTotal << endl;
    int idx = sprintf(buf, "%d\n", pathTotal);
    buf[idx] = '\0';
    fwrite(buf, idx, sizeof(char), fp);

    int outCountTotal = 0;
    int root = 0;
    int outPath3Count[4] = {0};
    int outPath4Count[4] = {0};
    int outPath5Count[4] = {0};
    int outPath6Count[4] = {0};
    int outPath7Count[4] = {0};
    int sos = 0;
    int outc = 0;

    root = 0;
    while (root < tnCount) {
      for (int i = 0; i < 4; i++) {
        while (path3[i][outPath3Count[i]][0] == root &&
               path3[i][outPath3Count[i]][1] != 0) {
          for (int j = 0; j < 2; j++) {
            memcpy(out + sos, idsComma[path3[i][outPath3Count[i]][j]].c_str(),
                   idsComma[path3[i][outPath3Count[i]][j]].size());
            sos += idsComma[path3[i][outPath3Count[i]][j]].size();
          }
          memcpy(out + sos, idsLF[path3[i][outPath3Count[i]][2]].c_str(),
                 idsLF[path3[i][outPath3Count[i]][2]].size());
          sos += idsLF[path3[i][outPath3Count[i]][2]].size();
          outc++;
          if (outc > 20000) {
            fwrite(out, sos, sizeof(char), fp);
            outc = 0;
            sos = 0;
          }
          outPath3Count[i]++;
        }
        root++;
      }
    }
    if (outc != 0) {
      fwrite(out, sos, sizeof(char), fp);
      outc = 0;
      sos = 0;
    }

    root = 0;
    while (root < tnCount) {
      for (int i = 0; i < 4; i++) {
        while (path4[i][outPath4Count[i]][0] == root &&
               path4[i][outPath4Count[i]][1] != 0) {
          for (int j = 0; j < 3; j++) {
            memcpy(out + sos, idsComma[path4[i][outPath4Count[i]][j]].c_str(),
                   idsComma[path4[i][outPath4Count[i]][j]].size());
            sos += idsComma[path4[i][outPath4Count[i]][j]].size();
          }
          memcpy(out + sos, idsLF[path4[i][outPath4Count[i]][3]].c_str(),
                 idsLF[path4[i][outPath4Count[i]][3]].size());
          sos += idsLF[path4[i][outPath4Count[i]][3]].size();
          outc++;
          if (outc > 15000) {
            fwrite(out, sos, sizeof(char), fp);
            outc = 0;
            sos = 0;
          }
          outPath4Count[i]++;
        }
        root++;
      }
    }
    if (outc != 0) {
      fwrite(out, sos, sizeof(char), fp);
      outc = 0;
      sos = 0;
    }

    root = 0;
    while (root < tnCount) {
      for (int i = 0; i < 4; i++) {
        while (path5[i][outPath5Count[i]][0] == root &&
               path5[i][outPath5Count[i]][1] != 0) {
          for (int j = 0; j < 4; j++) {
            memcpy(out + sos, idsComma[path5[i][outPath5Count[i]][j]].c_str(),
                   idsComma[path5[i][outPath5Count[i]][j]].size());
            sos += idsComma[path5[i][outPath5Count[i]][j]].size();
          }
          memcpy(out + sos, idsLF[path5[i][outPath5Count[i]][4]].c_str(),
                 idsLF[path5[i][outPath5Count[i]][4]].size());
          sos += idsLF[path5[i][outPath5Count[i]][4]].size();
          outc++;
          if (outc > 10000) {
            fwrite(out, sos, sizeof(char), fp);
            outc = 0;
            sos = 0;
          }
          outPath5Count[i]++;
        }
        root++;
      }
    }
    if (outc != 0) {
      fwrite(out, sos, sizeof(char), fp);
      outc = 0;
      sos = 0;
    }

    root = 0;
    while (root < tnCount) {
      for (int i = 0; i < 4; i++) {
        while (path6[i][outPath6Count[i]][0] == root &&
               path6[i][outPath6Count[i]][1] != 0) {
          for (int j = 0; j < 5; j++) {
            memcpy(out + sos, idsComma[path6[i][outPath6Count[i]][j]].c_str(),
                   idsComma[path6[i][outPath6Count[i]][j]].size());
            sos += idsComma[path6[i][outPath6Count[i]][j]].size();
          }
          memcpy(out + sos, idsLF[path6[i][outPath6Count[i]][5]].c_str(),
                 idsLF[path6[i][outPath6Count[i]][5]].size());
          sos += idsLF[path6[i][outPath6Count[i]][5]].size();
          outc++;
          if (outc > 10000) {
            fwrite(out, sos, sizeof(char), fp);
            outc = 0;
            sos = 0;
          }
          outPath6Count[i]++;
        }
        root++;
      }
    }
    if (outc != 0) {
      fwrite(out, sos, sizeof(char), fp);
      outc = 0;
      sos = 0;
    }

    root = 0;
    while (root < tnCount) {
      for (int i = 0; i < 4; i++) {
        while (path7[i][outPath7Count[i]][0] == root &&
               path7[i][outPath7Count[i]][1] != 0) {
          for (int j = 0; j < 6; j++) {
            memcpy(out + sos, idsComma[path7[i][outPath7Count[i]][j]].c_str(),
                   idsComma[path7[i][outPath7Count[i]][j]].size());
            sos += idsComma[path7[i][outPath7Count[i]][j]].size();
          }
          memcpy(out + sos, idsLF[path7[i][outPath7Count[i]][6]].c_str(),
                 idsLF[path7[i][outPath7Count[i]][6]].size());
          sos += idsLF[path7[i][outPath7Count[i]][6]].size();
          outc++;
          if (outc > 10000) {
            fwrite(out, sos, sizeof(char), fp);
            outc = 0;
            sos = 0;
          }
          outPath7Count[i]++;
        }
        root++;
      }
    }
    if (outc != 0) {
      fwrite(out, sos, sizeof(char), fp);
      outc = 0;
      sos = 0;
    }
    fclose(fp);
  }
};

int main() {
  string testFile = "./data/18908526/test_data.txt";
  string outputFile = "./data/18908526/result.txt";
  Solution solution;
  solution.parseInput(testFile);
  solution.constructGraph();
  solution.topologicalSort();
  solution.reConstructGraph();
  int t = clock();
  solution.ConstructPrecord();
  cout << "The cos time is:" << (double)(clock() - t) / CLOCKS_PER_SEC << "s"
       << endl;
  solution.solve();
  solution.save_fwrite(outputFile);
  cout << "The run time is:" << (double)clock() / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}