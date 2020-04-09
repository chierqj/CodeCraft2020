#include <arm_neon.h>
#include <fcntl.h>
#include <stdio.h>
#include <sys/ipc.h>
#include <sys/mman.h>
#include <sys/shm.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>
using namespace std;
using namespace std::chrono;

#define NOACCESS -1
#define HAVEACCESS -2
#define ACCESING -3
#define MMAX 280000

#define TEST
#define TIME

#ifdef TIME
steady_clock::time_point start = steady_clock::now();
#endif

#ifdef TEST
#define testFile "./data/test_data_10000_40000.txt"
#define resultFile "./data/result.txt"
#else
#define testFile "/data/test_data.txt"
#define resultFile "/projects/student/result.txt"
#endif
std::unordered_map<int, int> IDToMap;
int IDDom[MMAX];
int IDDomShadow[MMAX];

int IndexBack[MMAX];
int IndexForw[MMAX];

int Map[MMAX][3][100];
int childrenCount[MMAX];
int parentsCount[MMAX];

int circleTimes = 0;

int Circle[5][3500000][7];
int CircleCount[5];

struct dictOne {
  int start;
  char charNum[10];
};

void SaveAnswer(int maxID) {
  dictOne *dict = new dictOne[maxID];
  char *Answer = new char[circleTimes * 80 + 15];
  int AnswerIndex = 0;  // circleTimes*80+15-1;
  int Num = 0;
  for (int i = 0; i < maxID; i++) {
    Num = IDDom[i];
    if (Num == 0) {
      dict[i].charNum[9] = '0';
      dict[i].start = 9;
    } else {
      int j = 9;
      while (Num > 0) {
        dict[i].charNum[j--] = Num % 10 + '0';
        Num /= 10;
      }
      dict[i].start = j + 1;
    }
  }
  char CT[15];
  int CTNum = 14;
  Num = circleTimes;
  while (Num > 0) {
    CT[CTNum--] = Num % 10 + '0';
    Num /= 10;
  }
  CTNum++;
  for (int i = CTNum; i < 15; i++) {
    Answer[AnswerIndex++] = CT[i];
  }
  Answer[AnswerIndex++] = '\n';
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < CircleCount[i]; j++) {
      for (int k = 0; k < i + 3; k++) {
        for (int m = dict[Circle[i][j][k]].start; m < 10; m++) {
          Answer[AnswerIndex++] = dict[Circle[i][j][k]].charNum[m];
        }
        if (k != i + 2) {
          Answer[AnswerIndex++] = ',';
        }
      }
      Answer[AnswerIndex++] = '\n';
    }
  }
  Answer[AnswerIndex] = 0;
  FILE *fp = fopen(resultFile, "wt");
  fputs(Answer, fp);
  fclose(fp);
}

int LoadMap()  //维护一个邻接矩阵，一个父节点数目向量，一个子节点数目向量，返回值是最大ID
{
  int maxID = 0;
  for (int i = 0; i < MMAX; i++) {
    IndexForw[i] = NOACCESS;  //初始化为所有节点不存在
    childrenCount[i] = 0;     //子节点数目清空
    parentsCount[i] = 0;      //父节点数目清空
  }
  int oneRecord[3] = {0, 0, 0};
  char *p = NULL;
  char *data = NULL;
  int fd = open(testFile, O_RDONLY);
  long size = lseek(fd, 0, SEEK_END);
  data = (char *)mmap(NULL, size, PROT_READ, MAP_PRIVATE, fd, 0);
  p = data;
  int x = 0;
  while (p - data < size) {
    x = 0;
    while (*p != ',') {
      x *= 10;
      x += *p - '0';
      p++;
    }
    oneRecord[0] = x;
    x = 0;
    p++;
    while (*p != ',') {
      x *= 10;
      x += *p - '0';
      p++;
    }
    oneRecord[1] = x;
    x = 0;
    p++;
    while (*p != '\n') {
      if (*p >= 48 && *p <= 57) {
        x *= 10;
        x += *p - '0';
      }
      p++;
    }
    oneRecord[2] = x;
    p++;

    if (IDToMap.find(oneRecord[0]) == IDToMap.end()) {
      IDToMap[oneRecord[0]] = maxID;
      IDDom[maxID] = oneRecord[0];
      IDDomShadow[maxID] = oneRecord[0];
      maxID++;
    }
    if (IDToMap.find(oneRecord[1]) == IDToMap.end()) {
      IDToMap[oneRecord[1]] = maxID;
      IDDom[maxID] = oneRecord[1];
      IDDomShadow[maxID] = oneRecord[1];
      maxID++;
    }
    int one = IDToMap[oneRecord[1]];
    int zero = IDToMap[oneRecord[0]];
    Map[one][0][parentsCount[one]] = zero;
    Map[zero][1][childrenCount[zero]] = one;
    Map[zero][2][childrenCount[zero]] = oneRecord[2];
    childrenCount[zero]++;
    parentsCount[one]++;
  }
  return maxID;
}

void SortMap(int maxID) {
  for (int i = 0; i < maxID; i++) {
    if (childrenCount[i] <= 1)
      continue;
    else {
      for (int j = 0; j < childrenCount[i]; j++) {
        int minIndex = j;
        for (int k = minIndex + 1; k < childrenCount[i]; k++) {
          if (IDDom[Map[i][1][k]] < IDDom[Map[i][1][minIndex]]) {
            minIndex = k;
          }
        }
        int temp = Map[i][1][j];
        Map[i][1][j] = Map[i][1][minIndex];
        Map[i][1][minIndex] = temp;
        temp = Map[i][2][j];
        Map[i][2][j] = Map[i][2][minIndex];
        Map[i][2][minIndex] = temp;
      }
    }
  }
  for (int i = 0; i < maxID; i++) {
    if (parentsCount[i] <= 1)
      continue;
    else {
      for (int j = 0; j < parentsCount[i]; j++) {
        int minIndex = j;
        for (int k = minIndex + 1; k < parentsCount[i]; k++) {
          if (IDDom[Map[i][0][k]] < IDDom[Map[i][0][minIndex]]) {
            minIndex = k;
          }
        }
        int temp = Map[i][0][j];
        Map[i][0][j] = Map[i][0][minIndex];
        Map[i][0][minIndex] = temp;
      }
    }
  }
}

vector<pair<int, int>> BackPath[3];
int PathCount[3];
int NowRoot = 0;
void Back(int dot, int deep, int childern) {
  if (deep > 3 || IndexForw[dot] == HAVEACCESS || (dot == NowRoot && deep != 0))
    return;
  if (deep > 0) {
    BackPath[deep - 1].push_back(make_pair(dot, childern));
  }
  IndexBack[dot] = NowRoot;
  for (int i = 0; i < parentsCount[dot]; i++) {
    Back(Map[dot][0][i], deep + 1, dot);
  }
  return;
}
int chain[7];
bool sortPair(pair<int, int> &a, pair<int, int> &b) {
  if (a.first != b.first)
    return IDDom[a.first] < IDDom[b.first];
  else
    return IDDom[a.second] < IDDom[b.second];
}
bool sortInt(int &a, int &b) { return IDDom[a] < IDDom[b]; }
int ShortPath[500];
pair<int, int> LastPath[500];
int ShortCount = 0;
int LongCount = 0;
void Forward(int dot, int deep) {
  if (deep > 4 || IndexForw[dot] == HAVEACCESS || IndexForw[dot] == ACCESING ||
      (dot == NowRoot && deep != 0))
    return;
  IndexForw[dot] = ACCESING;
  if (IndexBack[dot] == NowRoot && deep > 0) {
    if (deep == 1) {
      ShortCount = 0;
      std::pair<int, int> faPair{dot, 0};
      std::pair<int, int> faPair1{dot, 560000};
      auto beg = std::lower_bound(BackPath[1].begin(),
                                  BackPath[1].begin() + PathCount[1], faPair);
      auto end = std::upper_bound(BackPath[1].begin(),
                                  BackPath[1].begin() + PathCount[1], faPair1);
      for (auto &it = beg; it != end; ++it) {
        ShortPath[ShortCount++] = it->second;
      }
      if (ShortCount > 0) {
        sort(ShortPath, ShortPath + ShortCount, sortInt);
        for (int i = 0; i < ShortCount; i++) {
          Circle[0][CircleCount[0]][0] = chain[0];
          Circle[0][CircleCount[0]][1] = dot;
          Circle[0][CircleCount[0]][2] = ShortPath[i];
          CircleCount[0]++;
          circleTimes++;
        }
      }
    }
    LongCount = 0;

    std::pair<int, int> faPair1{dot, 0};
    std::pair<int, int> faPair2{dot, 560000};
    auto beg1 = std::lower_bound(BackPath[2].begin(),
                                 BackPath[2].begin() + PathCount[2], faPair1);
    auto end1 = std::upper_bound(BackPath[2].begin(),
                                 BackPath[2].begin() + PathCount[2], faPair2);
    for (auto &it = beg1; it != end1; ++it) {
      std::pair<int, int> faPair3{it->second, -1};
      std::pair<int, int> faPair4{it->second, 560000};
      auto beg2 = std::lower_bound(BackPath[1].begin(),
                                   BackPath[1].begin() + PathCount[1], faPair3);
      auto end2 = std::upper_bound(BackPath[1].begin(),
                                   BackPath[1].begin() + PathCount[1], faPair4);
      for (auto &it2 = beg2; it2 != end2; ++it2) {
        if (IndexForw[it2->second] != ACCESING &&
            IndexForw[it2->first] != ACCESING) {
          LastPath[LongCount].first = it2->first;
          LastPath[LongCount].second = it2->second;
          LongCount++;
        }
      }
    }
    if (LongCount == 1) {
      int i = 0;
      for (i = 0; i < deep; i++) {
        Circle[deep][CircleCount[deep]][i] = chain[i];
      }
      Circle[deep][CircleCount[deep]][i++] = dot;
      Circle[deep][CircleCount[deep]][i++] = LastPath[0].first;
      Circle[deep][CircleCount[deep]][i] = LastPath[0].second;
      CircleCount[deep]++;
      circleTimes++;
    } else if (LongCount > 1) {
      sort(LastPath, LastPath + LongCount, sortPair);
      for (int j = 0; j < LongCount; j++) {
        int i = 0;
        for (i = 0; i < deep; i++) {
          Circle[deep][CircleCount[deep]][i] = chain[i];
        }
        Circle[deep][CircleCount[deep]][i++] = dot;
        Circle[deep][CircleCount[deep]][i++] = LastPath[j].first;
        Circle[deep][CircleCount[deep]][i] = LastPath[j].second;
        CircleCount[deep]++;
        circleTimes++;
      }
    }
  }
  for (int i = 0; i < childrenCount[dot]; i++) {
    chain[deep] = dot;
    Forward(Map[dot][1][i], deep + 1);
  }
  IndexForw[dot] = NOACCESS;
}

void BothWayFind(int maxID) {
  sort(IDDomShadow, IDDomShadow + maxID - 1);
  for (int i = 0; i < maxID; i++) {
    BackPath[0].clear(), BackPath[1].clear(), BackPath[2].clear();
    NowRoot = IDToMap[IDDomShadow[i]];
    Back(NowRoot, 0, NowRoot);
    for (int j = 0; j < 3; j++) {
      sort(BackPath[j].begin(), BackPath[j].end());
      PathCount[j] =
          unique(BackPath[j].begin(), BackPath[j].end()) - BackPath[j].begin();
    }
    Forward(NowRoot, 0);
    IndexForw[NowRoot] = HAVEACCESS;
  }
}

int main() {
#ifdef TIME
  cout << "开始加载邻接图.....";
  steady_clock::time_point startLoad = steady_clock::now();
#endif

  int maxID = LoadMap();
  maxID++;

#ifdef TIME
  cout << "节点数：" << maxID << ",时间："
       << duration_cast<duration<double>>(steady_clock::now() - startLoad)
              .count()
       << endl;
  cout << "开始预排序邻接图.....";
  steady_clock::time_point startSort = steady_clock::now();
#endif

  SortMap(maxID);

#ifdef TIME
  cout << "时间：";
  cout << duration_cast<duration<double>>(steady_clock::now() - startSort)
              .count()
       << endl;
  cout << "开始找环.....";
  steady_clock::time_point startFind = steady_clock::now();
#endif

  BothWayFind(maxID);

#ifdef TIME
  cout << "找到" << circleTimes << "个环,时间：";
  cout << duration_cast<duration<double>>(steady_clock::now() - startFind)
              .count()
       << endl;
  for (int i = 0; i < 5; i++) {
    cout << "长度" << i + 3 << "的环的数量是：" << CircleCount[i] << endl;
  }
  cout << "开始输出.....";
  steady_clock::time_point startOut = steady_clock::now();
#endif

  SaveAnswer(maxID);

#ifdef TIME
  cout << "输出结束,时间：";
  cout
      << duration_cast<duration<double>>(steady_clock::now() - startOut).count()
      << endl;
#endif
  return 0;
}
