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
#include <vector>

using namespace std;
using namespace std::chrono;
#define NOTEXIST -2
#define NOACCESS 0
#define ACCESSING 1
#define HAVEACCESS 2

#define TEST
#define TIME

#ifdef TIME
steady_clock::time_point start = steady_clock::now();
#endif

#ifdef TEST
#define testFile "../data/1004812/test_data.txt"
#define resultFile "../data/1004812/result.txt"
#else
#define testFile "/data/test_data_self.txt"
#define resultFile "/data/result.txt"
#endif
std::map<int, int> IDToMap;
std::set<int> AccessiblThree[280000];
std::set<int> AccessiblTwo[280000];
std::set<int> AccessiblOne[280000];
int IDDom[280000];
int IDDomShadow[280000];
int Index[280000];
int Map[280000][100][2];
int childrenCount[280000];
int parentsCount[280000];
int chain[7];
int circleTimes = 0;
int cutTimes = 0;

int AccessCount[280000];

int Circle[5][1000000][7];
int CircleCount[5];
char Answer[35000000];
int AnswerIndex = 34999999;

void SaveAnswer() {
  Answer[AnswerIndex--] = 0;
  for (int i = 4; i >= 0; i--) {
    for (int j = CircleCount[i] - 1; j >= 0; j--) {
      Answer[AnswerIndex--] = '\n';
      for (int k = i + 2; k >= 0; k--) {
        int Num = IDDom[Circle[i][j][k]];
        while (Num > 0) {
          Answer[AnswerIndex--] = Num % 10 + '0';
          Num /= 10;
        }
        if (k != 0) {
          Answer[AnswerIndex--] = ',';
        }
      }
    }
  }
  Answer[AnswerIndex--] = '\n';
  int Num = circleTimes;
  while (Num > 0) {
    Answer[AnswerIndex--] = Num % 10 + '0';
    Num /= 10;
  }
  char *answerPtr = &Answer[AnswerIndex + 1];
  FILE *fp = fopen(resultFile, "wt");
  fputs(answerPtr, fp);
  fclose(fp);
}

void Gao4(int dot, int deep) {
  // AccessCount[dot]++;
  if (deep > 6 || Index[dot] != NOACCESS) return;
  int root = chain[0];
  if (deep == 4) {
    if (AccessiblThree[dot].find(root) == AccessiblThree[dot].end()) return;
  } else if (deep == 5) {
    if (AccessiblTwo[dot].find(root) == AccessiblTwo[dot].end()) return;
  } else if (deep == 6) {
    if (AccessiblOne[dot].find(root) == AccessiblOne[dot].end()) return;
  }

  chain[deep] = dot;
  Index[dot] = ACCESSING;
  int Count = childrenCount[dot];
  for (int i = 0; i < Count; i++) {
    if (deep > 1 && Map[dot][i][0] == root) {
      for (int j = 0; j <= deep; j++) {
        Circle[deep + 1 - 3][CircleCount[deep + 1 - 3]][j] = chain[j];
      }
      CircleCount[deep + 1 - 3]++;
      circleTimes++;
      continue;
    } else {
      Gao4(Map[dot][i][0], deep + 1);
    }
  }
  Index[dot] = NOACCESS;
  return;
}

void Gao3(int dot, int farther, int deep) {
  if (deep == 0) {
    for (int i = 0; i < childrenCount[dot]; i++) {
      Gao3(Map[dot][i][0], farther, deep + 1);
    }
  } else if (deep == 1) {
    AccessiblOne[farther].insert(dot);
    AccessiblTwo[farther].insert(dot);
    AccessiblThree[farther].insert(dot);
    for (int i = 0; i < childrenCount[dot]; i++) {
      Gao3(Map[dot][i][0], farther, deep + 1);
    }
  } else if (deep == 2) {
    AccessiblTwo[farther].insert(dot);
    AccessiblThree[farther].insert(dot);
    for (int i = 0; i < childrenCount[dot]; i++) {
      Gao3(Map[dot][i][0], farther, deep + 1);
    }
  } else if (deep == 3) {
    AccessiblThree[farther].insert(dot);
    return;
  }

  return;
}

int LoadMap()  //维护一个邻接矩阵，一个父节点数目向量，一个子节点数目向量，返回值是最大ID
{
  int maxID = 0;
  int oneRecord[3] = {0, 0, 0};
  char c;
  char *p;
  char *data = NULL;
  int fd = open(testFile, O_RDONLY);
  long size = lseek(fd, 0, SEEK_END);
  data = (char *)mmap(NULL, size, PROT_READ, MAP_PRIVATE, fd, 0);
  p = data;
  while (p - data < size) {
    int x = 0;
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
    int count = childrenCount[zero];
    Map[zero][count][0] = one;
    Map[zero][count][1] = oneRecord[2];
    childrenCount[zero]++;
  }
  return maxID;
}

void SortMap(int maxID) {
  for (int i = 0; i < maxID; i++) {
    if (Index[i] == NOTEXIST || childrenCount[i] <= 1)
      continue;
    else {
      for (int j = 0; j < childrenCount[i]; j++) {
        int minIndex = j;
        for (int k = minIndex + 1; k < childrenCount[i]; k++) {
          if (IDDom[Map[i][k][1]] < IDDom[Map[i][minIndex][1]]) {
            minIndex = k;
          }
        }
        int temp = Map[i][j][1];
        Map[i][j][1] = Map[i][minIndex][1];
        Map[i][minIndex][1] = temp;
        temp = Map[i][j][2];
        Map[i][j][2] = Map[i][minIndex][2];
        Map[i][minIndex][2] = temp;
      }
    }
  }
  /*for(int i = 0;i<maxID;i++)
  {
          for(int j = 0;j<childrenCount[i];j++)
          {
                  cout<<IDDom[Map[i][j][1]]<<",";
          }
          cout<<endl;
  }*/
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
  cout << "Start SortMap: ";
  cout << duration_cast<duration<double>>(steady_clock::now() - start).count()
       << endl;
#endif

  SortMap(maxID);

#ifdef TIME
  cout << "Start DeFind: ";
  cout << duration_cast<duration<double>>(steady_clock::now() - start).count()
       << endl;
#endif

  for (int i = 0; i < maxID; i++) {
    if (Index[i] == NOACCESS) {
      Gao3(i, i, 0);
    }
  }
#ifdef TIME
  cout << "Start Find: ";
  cout << duration_cast<duration<double>>(steady_clock::now() - start).count()
       << endl;
#endif
  sort(IDDomShadow, IDDomShadow + maxID - 1);
  for (int i = 0; i < maxID; i++) {
    if (Index[IDToMap[IDDomShadow[i]]] == NOACCESS) {
      Gao4(IDToMap[IDDomShadow[i]], 0);
      Index[IDToMap[IDDomShadow[i]]] = HAVEACCESS;
    }
  }
  int MeanCount = 0, MeanParents = 0;
  ;
  for (int i = 0; i < maxID; i++) {
    MeanCount += AccessCount[i];
    MeanParents += parentsCount[i];
    // cout<<i<<"AccessCount: "<<AccessCount[i]<<".Prents:
    // "<<parentsCount[i]<<endl;
  }
  cout << "MeanAccessCount: " << MeanCount / maxID
       << ".MeanParents: " << MeanParents / maxID << endl;
#ifdef TIME
  cout << "\t\t\t\tCircle: " << circleTimes << endl;
  cout << "Start Out: ";
  cout << duration_cast<duration<double>>(steady_clock::now() - start).count()
       << endl;
#endif
  SaveAnswer();

#ifdef TIME
  cout << "All Time: ";
  cout << duration_cast<duration<double>>(steady_clock::now() - start).count()
       << endl;
#endif
  return 0;
}