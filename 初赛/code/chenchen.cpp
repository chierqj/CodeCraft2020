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
#include <memory>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace std;
using namespace std::chrono;
#define NOTEXIST -2
#define NOACCESS -1
#define ACCESSING 1
#define HAVEACCESS 2

#define TEST
#define TIME

#ifdef TIME
steady_clock::time_point start = steady_clock::now();
#endif

#ifdef TEST
#define testFile "./data/test_data2.txt"
#define resultFile "./data/result.txt"
#else
#define testFile "/data/test_data_self.txt"
#define resultFile "/data/result.txt"
#endif
int Index[280000];
int Map[280000][100][3];
int childrenCount[280000];
int parentsCount[280000];
int chain[15000];
vector<vector<vector<int>>> Circle;
int circleTimes = 0;
int cutTimes = 0;
void Gao(int dot, int indexNum) {
  if (dot == NOTEXIST || Index[dot] == NOTEXIST || Index[dot] == HAVEACCESS ||
      indexNum >= 7 ||
      Index[dot] == ACCESSING) {  //判断是否存在，是否已经找过环
    return;
  } else {
    chain[indexNum] = dot;
    Index[dot] = ACCESSING;
    for (int i = 0; i < childrenCount[dot]; i++)  //遍历子节点
    {
      if (Map[dot][i][1] == chain[0])  //如果子节点是父节点
      {
        if (indexNum + 1 >= 3 && indexNum + 1 <= 7)  //如果长度合适
        {
          vector<int> Record;
          for (int j = 0; j <= indexNum; j++) {
            Record.push_back(chain[j]);
          }
          Circle[indexNum + 1 - 3].push_back(Record);
          circleTimes++;
        }

      } else {
        Gao(Map[dot][i][1], indexNum + 1);
      }
    }
    Index[dot] = NOACCESS;
  }
  return;
}

int LoadMap()  //维护一个邻接矩阵，一个父节点数目向量，一个子节点数目向量，返回值是最大ID
{
  int maxID = 0;
  for (int i = 0; i < 280000; i++) {
    Index[i] = NOTEXIST;   //初始化为所有节点不存在
    childrenCount[i] = 0;  //子节点数目清空
    parentsCount[i] = 0;   //父节点数目清空
  }
  int oneRecord[3] = {0, 0, 0};
  char c;
  int8_t *p;
  int8_t *data = NULL;
  int fd = open(testFile, O_RDONLY);
  long size = lseek(fd, 0, SEEK_END);
  data = (int8_t *)mmap(NULL, size, PROT_READ, MAP_PRIVATE, fd, 0);
  p = data;
  while (true) {
    while (*p != ',') {
      oneRecord[0] *= 10;
      oneRecord[0] += *p - '0';
      p++;
    }
    p++;
    while (*p != ',') {
      oneRecord[1] *= 10;
      oneRecord[1] += *p - '0';
      p++;
    }
    p++;
    while (*p != '\n') {
      if (*p >= 48 && *p <= 57) {
        oneRecord[2] *= 10;
        oneRecord[2] += *p - '0';
      }
      p++;
    }
    p++;
    if (Index[oneRecord[0]] == NOTEXIST) {
      Index[oneRecord[0]] = NOACCESS;
    }
    if (Index[oneRecord[1]] == NOTEXIST) {
      Index[oneRecord[1]] = NOACCESS;
    }
    if (oneRecord[0] > maxID)
      maxID = oneRecord[0];
    else if (oneRecord[1] > maxID)
      maxID = oneRecord[1];
    Map[oneRecord[1]][parentsCount[oneRecord[1]]][0] = oneRecord[0];
    Map[oneRecord[0]][childrenCount[oneRecord[0]]][1] = oneRecord[1];
    Map[oneRecord[0]][childrenCount[oneRecord[0]]][2] = oneRecord[2];
    childrenCount[oneRecord[0]]++;
    parentsCount[oneRecord[1]]++;
    oneRecord[0] = 0;
    oneRecord[1] = 0;
    oneRecord[2] = 0;
    if (p - data >= size) break;
  }
  return maxID;
}

void InitCircle(void)  //初始化保存环的向量，不同长度的环存的位置不同
{
  vector<int> circle1;
  vector<vector<int>> circleHome;
  circle1.push_back(-1);
  circleHome.push_back(circle1);
  for (int i = 0; i < 5; i++) {
    Circle.push_back(circleHome);
  }
}

void forwardRecursion(int dot, int parents) {
  if (parents != NOTEXIST && parentsCount[dot] != 0) {
    for (int i = 0; i < parentsCount[dot]; i++) {
      if (Map[dot][i][0] == parents) Map[dot][i][0] = NOTEXIST;
    }
    parentsCount[dot]--;
  }

  if (childrenCount[dot] != 0 && parentsCount[dot] == 0)  //被剪枝后成根节点
  {
    for (int i = 0; i < childrenCount[dot]; i++) {
      forwardRecursion(Map[dot][i][1], dot);
    }
    Index[dot] = NOTEXIST;
    cutTimes++;
  } else if (childrenCount[dot] == 0 &&
             parentsCount[dot] == 0)  //被剪枝后成孤立节点
  {
    Index[dot] = NOTEXIST;
    cutTimes++;
  }
  return;
}

void backwardRecursion(int dot, int children) {
  if (children != NOTEXIST && childrenCount[dot] != 0) {
    for (int i = 0; i < childrenCount[dot]; i++) {
      if (Map[dot][i][1] == children) Map[dot][i][1] = NOTEXIST;
    }
    childrenCount[dot]--;
  }

  if (childrenCount[dot] == 0 && parentsCount[dot] != 0)  //被剪枝后成根节点
  {
    for (int i = 0; i < parentsCount[dot]; i++) {
      forwardRecursion(Map[dot][i][0], dot);
    }
    Index[dot] = NOTEXIST;
    cutTimes++;
  } else if (childrenCount[dot] == 0 &&
             parentsCount[dot] == 0)  //被剪枝后成孤立节点
  {
    Index[dot] = NOTEXIST;
    cutTimes++;
  }
  return;
}

void cutCuteBranch(int maxID) {
  for (int i = 0; i < maxID; i++) {
    if (Index[i] != NOTEXIST) {
      forwardRecursion(i, NOTEXIST);  //子节点是i,不输入父节点
      backwardRecursion(i, NOTEXIST);
    }
  }
}

int main() {
#ifdef TIME
  cout << "Start LoadMap: ";
  cout << duration_cast<duration<double>>(steady_clock::now() - start).count()
       << endl;
#endif
  int maxID = LoadMap();
  maxID++;
#ifdef TIME
  cout << "Start CutBranch: ";
  cout << duration_cast<duration<double>>(steady_clock::now() - start).count()
       << endl;
#endif

  cutCuteBranch(maxID);
  cout << "\t\t\t\tCutBranch: " << cutTimes << endl;
  InitCircle();

#ifdef TIME
  cout << "Start Find: ";
  cout << duration_cast<duration<double>>(steady_clock::now() - start).count()
       << endl;
#endif

  for (int i = 0; i <= maxID; i++) {
    if (Index[i] == NOACCESS) {
      Gao(i, 0);
      Index[i] = HAVEACCESS;
      cout << i << " " << circleTimes << endl;
    }
  }

#ifdef TIME
  cout << "\t\t\t\tCircle: " << circleTimes << endl;
  cout << "Start Sort: ";
  cout << duration_cast<duration<double>>(steady_clock::now() - start).count()
       << endl;
#endif
  for (int i = 0; i < 5; i++) {
    std::sort(Circle[i].begin(), Circle[i].end(),
              [&](const std::vector<int> &vt1, const std::vector<int> &vt2) {
                return vt1 < vt2;
              });
  }
/*#ifdef TIME
cout << "Start Out: ";
cout << duration_cast<duration<double>>(steady_clock::now() - start).count() <<
endl; #endif ofstream fout(resultFile); fout<<circleTimes<<endl; for(int i =
0;i<5;i++)
{
        for(int j = 1;j<Circle[i].size();j++)
        {
                fout<<Circle[i][j][0];
                for(int k = 1;k<Circle[i][j].size();k++)
                {
                        fout<<","<<Circle[i][j][k];
                }
                fout<<endl;
        }
}
fout.close(); */
#ifdef TIME
  cout << "All Time: ";
  cout << duration_cast<duration<double>>(steady_clock::now() - start).count()
       << endl;
#endif
  return 0;
}
