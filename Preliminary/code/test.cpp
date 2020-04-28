#include <bits/stdc++.h>
#include <time.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>  //lock
#include <numeric>
#include <queue>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>
using namespace std;
/*
±àÒë·½Ê½£ºg++ -std=c++11 -O3 Version23.cpp -o Version23 -lpthread
±àÒë·½Ê½£ºg++ -O3 Version23.cpp -o Version23 -lpthread
*/
struct Path {
  vector<int> path;
  Path(const vector<int> &path) : path(path) {}
};

#define DEPTH_LIMIT_HIGH 7
#define DEPTH_LIMIT_LOW 3
#define ROW 280000
#define COLUMN 50
class Solution {
 public:
  int count = 0;
  mutex m_mutex;
  vector<int> inputs;
  int nodeCnt;
  vector<int> ids;
  int **GG;

  unordered_map<int, int> idHash;  // sorted id to 0...n
  vector<int> inDegrees;
  vector<int> outDegrees;
  vector<unordered_map<int, vector<int>>> ring_three;
  vector<int> single_result_one;
  // vector<int> single_result_two;
  // vector<int> single_result_three;
  // vector<int> single_result_four;
  vector<vector<int>> result_one;
  // vector<vector<int>> result_two;
  // vector<vector<int>> result_three;
  // vector<vector<int>> result_four;
  vector<int> reachable_one;
  // vector<int> reachable_two;
  // vector<int> reachable_three;
  // vector<int> reachable_four;
  vector<int> cnode_one;
  // vector<int> cnode_two;
  // vector<int> cnode_three;
  // vector<int> cnode_four;
  // vector<bool> visited_one;
  bool *visited_one;
  // vector<bool> visited_two;
  // vector<bool> visited_three;
  // vector<bool> visited_four;
  int depth_one = 1;
  // int depth_two = 1;
  // int depth_three = 1;
  // int depth_four = 1;

  vector<string> idsComma;  // 0...n to sorted id
  vector<string> idsLF;     // 0...n to sorted id
  vector<Path> ans[8];

 public:
  Solution() : m_mutex() {}
  void parseInput(string &testFile) {
    FILE *file = fopen(testFile.c_str(), "r");
    int u, v, c;
    while (fscanf(file, "%d,%d,%d", &u, &v, &c) != EOF) {
      inputs.emplace_back(u);
      inputs.emplace_back(v);
    }
  }
  void constructGraph() {
    auto tmp = inputs;
    sort(tmp.begin(), tmp.end());
    tmp.erase(unique(tmp.begin(), tmp.end()), tmp.end());
    nodeCnt = tmp.size();
    idsComma.reserve(nodeCnt);
    idsLF.reserve(nodeCnt);
    nodeCnt = 0;
    for (int &x : tmp) {
      idsComma.push_back(to_string(x) + ',');
      idsLF.push_back(to_string(x) + '\n');
      idHash[x] = nodeCnt++;
    }
    int sz = inputs.size();
    GG = new int *[nodeCnt];
    for (int i = 0; i < nodeCnt; i++) {
      GG[i] = new int[COLUMN];
      memset(GG[i], -1, COLUMN * sizeof(int));
    }
    inDegrees = vector<int>(nodeCnt, 0);
    outDegrees = vector<int>(nodeCnt, 0);
    for (int i = 0; i < sz; i += 2) {
      int u = idHash[inputs[i]], v = idHash[inputs[i + 1]];
      GG[u][outDegrees[u]] = v;
      ++inDegrees[v];
      ++outDegrees[u];
    }
  }

  void constructGraphTwo() {
    ring_three = vector<unordered_map<int, vector<int>>>(
        nodeCnt, unordered_map<int, vector<int>>());
    for (int i = 0; i < nodeCnt; i++) {
      auto &one = GG[i];
      for (int j = 0; j < COLUMN && one[j] != -1; j++) {
        auto &two = GG[one[j]];
        for (int k = 0; k < COLUMN && two[k] != -1; k++) {
          if (two[k] != i) {
            ring_three[two[k]][i].push_back(one[j]);
          }
        }
      }
    }

    for (int i = 0; i < nodeCnt; i++) {
      for (auto &x : ring_three[i]) {
        if (x.second.size() > 1) {
          sort(x.second.begin(), x.second.end());
        }
      }
    }
  }
  void topoSort(vector<int> &degs, bool doSoring) {
    queue<int> que;
    for (int i = 0; i < nodeCnt; i++) {
      if (0 == degs[i]) que.push(i);
    }
    while (!que.empty()) {
      int u = que.front();
      que.pop();
      auto one = GG[u];
      for (int i = 0; i < COLUMN && one[i] != -1; i++) {
        auto v = one[i];
        if (0 == --degs[v]) {
          que.push(v);
        }
      }
    }
    int cnt = 0;
    for (int i = 0; i < nodeCnt; i++) {
      if (degs[i] == 0) {
        memset(GG[i], -1, sizeof(GG[i]));
        cnt++;
      } else if (doSoring) {
        sort(GG[i], GG[i] + outDegrees[i]);
      }
    }
  }

  void IterationAlgorithm(int node, bool *&visited, vector<int> &single_result,
                          int current_node, int depth) {
    visited[node] = true;
    single_result.push_back(node);
    auto gCur = GG[node];
    auto it = lower_bound(gCur, gCur + outDegrees[node], current_node);
    // auto it = lower_bound(gCur.begin(), gCur.end(), current_node);
    if (*it == current_node && depth >= DEPTH_LIMIT_LOW &&
        depth < DEPTH_LIMIT_HIGH) {
      // m_mutex.lock();
      this->count++;
      // m_mutex.unlock();
      ans[depth].emplace_back(Path(single_result));
    }

    if (depth < DEPTH_LIMIT_HIGH - 1) {
      for (; *it != -1; ++it) {
        if (!visited[*it]) {
          IterationAlgorithm(*it, visited, single_result, current_node,
                             depth + 1);
        }
      }
    } else if (reachable_one[node] > -1 && depth == DEPTH_LIMIT_HIGH - 1) {
      auto ks = ring_three[current_node][node];
      int sz = ks.size();
      for (int j = reachable_one[node]; j < sz; j++) {
        int k = ks[j];
        if (visited[k]) continue;
        single_result.push_back(k);
        // m_mutex.lock();
        this->count++;
        // m_mutex.unlock();
        // result_one.emplace_back(single_result);
        ans[depth + 1].emplace_back(Path(single_result));
        single_result.pop_back();
      }
    }
    visited[node] = false;
    single_result.pop_back();
  }
  void seven_reachable_one(int node) {
    for (auto &js : ring_three[node]) {
      int j = js.first;
      if (j > node) {
        auto &val = js.second;
        int lb = lower_bound(val.begin(), val.end(), node) - val.begin();
        if (lb < val.size()) reachable_one[j] = lb;
        cnode_one.push_back(j);
      }
    }
  }
  void seven_reachable_two() {
    for (int &x : cnode_one) reachable_one[x] = -1;
    cnode_one.clear();
  }
  void user_thread(int start, int end, int depth, bool *&visited,
                   vector<int> &single_result) {
    for (int i = start; i < end; i++) {
      seven_reachable_one(i);
      single_result.emplace_back(i);
      visited[i] = true;
      auto &one = GG[i];

      for (int j = 0; j < COLUMN && one[j] != -1; j++) {
        if (!visited[one[j]] && one[j] > i) {
          visited[one[j]] = true;
          depth++;
          single_result.emplace_back(one[j]);
          auto &two = GG[one[j]];
          for (int k = 0; k < COLUMN && two[k] != -1; k++) {
            if (!visited[two[k]] && two[k] > i) {
              depth++;
              IterationAlgorithm(two[k], ref(visited), ref(single_result), i,
                                 ref(depth));
              depth--;
            }
          }
          visited[one[j]] = false;
          single_result.pop_back();
          depth--;
        }
      }
      visited[i] = false;
      single_result.pop_back();
      seven_reachable_two();
    }
  }
  // void merge_thread(vector<vector<int>> &final_sresult){
  //	int size = final_sresult.size();
  //	for (int i = 0; i <size; i++){
  //		ans[final_sresult[i].size()].emplace_back(Path(final_sresult[i]));
  //	}
  //}
  void solve() {
    visited_one = new bool[nodeCnt]();
    // visited_two = vector<bool>(nodeCnt, false);
    // visited_three = vector<bool>(nodeCnt, false);
    // visited_four = vector<bool>(nodeCnt, false);

    reachable_one = vector<int>(nodeCnt, -1);
    // reachable_two = vector<int>(nodeCnt, -1);
    // reachable_three = vector<int>(nodeCnt, -1);
    // reachable_four = vector<int>(nodeCnt, -1);
    cnode_one = vector<int>(nodeCnt);
    // cnode_two = vector<int>(nodeCnt);
    // cnode_three = vector<int>(nodeCnt);
    // cnode_four = vector<int>(nodeCnt);
    // int node_number = nodeCnt;
    // int average_number = (node_number / 4) + 1;
    // thread task01(&Solution::user_thread, this, 0, nodeCnt, ref(depth_one),
    // ref(visited_one), 	ref(single_result_one), ref(reachable_one),
    //ref(cnode_one), ref(result_one));
    user_thread(0, nodeCnt, ref(depth_one), ref(visited_one),
                ref(single_result_one));
    // thread task02(&Solution::user_thread, this, average_number, 2 *
    // average_number, ref(depth_two), ref(visited_two), 	ref(single_result_two),
    //ref(reachable_two), ref(cnode_two), ref(result_two)); thread
    // task03(&Solution::user_thread, this, 2 * average_number, 3 *
    // average_number, ref(depth_three), ref(visited_three),
    //	ref(single_result_three), ref(reachable_three), ref(cnode_three),
    //ref(result_three)); thread task04(&Solution::user_thread, this, 3 *
    // average_number, nodeCnt, ref(depth_four), ref(visited_four),
    //	ref(single_result_four), ref(reachable_four), ref(cnode_four),
    //ref(result_four)); task01.join(); merge_thread(result_one); task02.join();
    // merge_thread(result_two);
    // task03.join();
    // merge_thread(result_three);
    // task04.join();
    // merge_thread(result_four);
  }
  void save_fwrite(const string &outputFile) {
    FILE *fp = fopen(outputFile.c_str(), "wb");
    char buf[1024];
    int idx = sprintf(buf, "%d\n", count);
    buf[idx] = '\0';
    fwrite(buf, idx, sizeof(char), fp);
    for (int i = DEPTH_LIMIT_LOW; i <= DEPTH_LIMIT_HIGH; i++) {
      for (auto &x : ans[i]) {
        auto path = x.path;
        int sz = path.size();
        for (int j = 0; j < sz - 1; j++) {
          auto res = idsComma[path[j]];
          fwrite(res.c_str(), res.size(), sizeof(char), fp);
        }
        auto res = idsLF[path[sz - 1]];
        fwrite(res.c_str(), res.size(), sizeof(char), fp);
      }
    }
    fclose(fp);
  }
};

int main() {
  // clock_t start, end;
  // start = clock();
  // string testFile = "./test_data3.txt";
  // string outputFile = "./output.txt";
  string testFile = "/data/test_data.txt";
  string outputFile = "/projects/student/result.txt";
  Solution sol;
  sol.parseInput(testFile);
  sol.constructGraph();
  sol.constructGraphTwo();
  sol.topoSort(sol.inDegrees, true);
  // sol.topoSort(sol.outDegrees, true);
  sol.solve();
  sol.save_fwrite(outputFile);
  // end = clock();
  // cout << "Time is " << (float)(end - start) / CLOCKS_PER_SEC << endl;
  return 0;
}