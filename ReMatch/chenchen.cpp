#include <arm_neon.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
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

#define MMAX 4000000
#define CHILDRENCORE 3  //子进程数目

unsigned int ChildrenTimes[MMAX];  //子节点数目
unsigned int ParentsTimes[MMAX];   //父节点数目
bool IndexBack[MMAX][3];
bool IndexForw[MMAX];
unsigned int Buff[MMAX][3];
unsigned int IDTable[MMAX * 2];  //存放预读取的ID
unsigned int circleTimes = 0;
unsigned int NowBuff = 0;
#define TEST
#define TIME
unsigned int ChildrenNum[MMAX];
unsigned int ParentNum[MMAX];
pair<unsigned int, unsigned int>* ChildrenTable;
pair<unsigned int, unsigned int>* ChildrenStart[MMAX];
char* CircleAnswer[5];
unsigned int CircleAnswerID[9][2];
char* CircleAnswerPtr[9][5];  //存的是9个开始点
long CircleAnswerIdx[9][5];   //存的是9个长度
unsigned int* ParentTable;
unsigned int* ParentStart[MMAX];
#ifdef TIME
steady_clock::time_point start = steady_clock::now();
#endif

#ifdef TEST
#define testFile "./data/19630345/test_data.txt"
#define resultFile "./data/19630345/result.txt"
#else
#define testFile "/data/test_data.txt"
#define resultFile "/projects/student/result.txt"
#endif
struct dictOne {
  unsigned int start;
  unsigned int length;
  char charNum[10];
};
dictOne* dict;
void BuildDict(unsigned int IDNum) {
  for (unsigned int i = 0; i < IDNum; i++)  //字典自己处理自己的
  {
    unsigned int Num = IDTable[i];
    if (Num == 0) {
      dict[i].charNum[9] = '0';
      dict[i].start = 9;
      dict[i].length = 1;
    } else {
      unsigned int j = 9;
      while (Num > 0) {
        dict[i].charNum[j--] = Num % 10 + '0';
        Num /= 10;
      }
      dict[i].start = j + 1;
      dict[i].length = 10 - dict[i].start;
    }
  }
  return;
}

void SaveAnswer(unsigned int maxID, unsigned int AllCircleTimes,
                unsigned int ChildrenNum, unsigned int* OUTPtr) {
  char* Answer = new char[15];
  unsigned int AnswerIndex = 0;
  unsigned int Num = 0;

  if (ChildrenNum == 0)  //由父节点写入环的数目
  {
    char CT[15];
    unsigned int CTNum = 14;
    Num = AllCircleTimes;
    while (Num > 0) {
      CT[CTNum--] = Num % 10 + '0';
      Num /= 10;
    }
    CTNum++;
    for (unsigned int i = CTNum; i < 15; i++) {
      Answer[AnswerIndex++] = CT[i];
    }
    Answer[AnswerIndex++] = '\n';
    Answer[AnswerIndex] = 0;
    FILE* fp = fopen(resultFile, "wt");
    fputs(Answer, fp);
    fclose(fp);
    *OUTPtr = 0;
  }
  for (unsigned int i = 0; i < 5; i++)  //输出五种环
  {
    for (unsigned int j = 0; j < NowBuff; j++) {
      while (true) {
        if (*OUTPtr == CircleAnswerID[j][0]) {
          break;
        }
        usleep(1);
      }
      FILE* fp = fopen(resultFile, "at");
      fwrite(CircleAnswerPtr[j][i], 1, CircleAnswerIdx[j][i], fp);
      fclose(fp);
      *OUTPtr = CircleAnswerID[j][1];
      if (*OUTPtr >= maxID) *OUTPtr = 0;
    }
  }
}

unsigned int
LoadMap()  //维护一个邻接矩阵，一个父节点数目向量，一个子节点数目向量，返回值是最大ID
{
  for (unsigned int i = 0; i < MMAX; i++) {
    ChildrenTimes[i] = 0;
    ParentsTimes[i] = 0;
    ChildrenNum[i] = 0;
    ParentNum[i] = 0;
  }
  for (unsigned int i = 0; i < 5; i++) {
    CircleAnswer[i] = new char[20000000 * (i + 3) * 11];
    for (unsigned int j = 0; j < 9; j++) {
      CircleAnswerIdx[j][i] = 0;
    }
  }
  char* p = NULL;
  char* data = NULL;
  unsigned int fd = open(testFile, O_RDONLY);
  long size = lseek(fd, 0, SEEK_END);
  data = (char*)mmap(NULL, size, PROT_READ, MAP_PRIVATE, fd, 0);
  p = data;
  unsigned int x = 0;
  unsigned int Index = 0;
  unsigned int IDIndex = 0;
  while (p - data < size) {
    x = 0;
    while (*p != ',') {
      x *= 10;
      x += *p - '0';
      p++;
    }
    Buff[Index][0] = x;
    IDTable[IDIndex++] = x;
    x = 0;
    p++;
    while (*p != ',') {
      x *= 10;
      x += *p - '0';
      p++;
    }
    Buff[Index][1] = x;
    IDTable[IDIndex++] = x;
    x = 0;
    p++;
    while (*p != '\n') {
      if (*p >= 48 && *p <= 57) {
        x *= 10;
        x += *p - '0';
      }
      p++;
    }
    Buff[Index][2] = x;
    p++;
    Index++;
  }
  sort(IDTable, IDTable + IDIndex);
  unsigned int IDNum = unique(IDTable, IDTable + IDIndex) - IDTable;
  /*for(unsigned int i = 0;i<20;i++)
  {
          cout<<IDTable[i]<<endl;
  }*/
  dict = new dictOne[IDNum];
  BuildDict(IDNum);
  ChildrenTable = new pair<unsigned int, unsigned int>[Index];
  ParentTable = new unsigned int[Index];
  for (unsigned int i = 0; i < Index; i++) {
    Buff[i][0] = lower_bound(IDTable, IDTable + IDNum, Buff[i][0]) - IDTable;
    Buff[i][1] = lower_bound(IDTable, IDTable + IDNum, Buff[i][1]) - IDTable;
    ChildrenTimes[Buff[i][0]]++;
    ParentsTimes[Buff[i][1]]++;
  }
  ChildrenStart[0] = ChildrenTable;
  ParentStart[0] = ParentTable;
  for (unsigned int i = 1; i < IDNum; i++) {
    ChildrenStart[i] = ChildrenStart[i - 1] + ChildrenTimes[i - 1];
    ParentStart[i] = ParentStart[i - 1] + ParentsTimes[i - 1];
  }
  for (unsigned int i = 0; i < Index; i++) {
    (ChildrenStart[Buff[i][0]] + ChildrenNum[Buff[i][0]])->first = Buff[i][1];
    (ChildrenStart[Buff[i][0]] + ChildrenNum[Buff[i][0]])->second = Buff[i][2];
    *(ParentStart[Buff[i][1]] + ParentNum[Buff[i][1]]) = Buff[i][0];
    ChildrenNum[Buff[i][0]] += 1;
    ParentNum[Buff[i][1]] += 1;
  }
  for (unsigned int i = 0; i < IDNum; i++) {
    sort(ChildrenStart[i], ChildrenStart[i] + ChildrenNum[i]);
    sort(ParentStart[i], ParentStart[i] + ParentNum[i]);
  }
  /*cout<<endl;
  for(unsigned int i = 0;i<ChildrenNum[0];i++)
  {
          cout<<(ChildrenStart[0]+i)->first<<endl;
  }*/
  return IDNum;
}
unsigned int Flag[MMAX];
unsigned int Count = 0;

void Back(unsigned int NowRoot) {
  unsigned int ccount1 = ParentNum[NowRoot];
  for (unsigned int i = 0; i < ccount1; i++) {
    unsigned int deep1dot = *(ParentStart[NowRoot] + i);  // Map[NowRoot][0][i];
    if (IndexForw[deep1dot]) continue;
    Flag[Count++] = deep1dot;
    IndexBack[deep1dot][0] = true;
    IndexBack[deep1dot][1] = true;
    IndexBack[deep1dot][2] = true;
    unsigned int ccount2 = ParentNum[deep1dot];
    for (unsigned int j = 0; j < ccount2; j++) {
      unsigned int deep2dot =
          *(ParentStart[deep1dot] + j);  // Map[deep1dot][0][j];
      if (IndexForw[deep2dot]) continue;
      Flag[Count++] = deep2dot;
      IndexBack[deep2dot][1] = true;
      IndexBack[deep2dot][2] = true;
      unsigned int ccount3 = ParentNum[deep2dot];
      for (unsigned int k = 0; k < ccount3; k++) {
        unsigned int deep3dot =
            *(ParentStart[deep2dot] + k);  // Map[deep2dot][0][k];
        if (IndexForw[deep3dot]) continue;
        Flag[Count++] = deep3dot;
        IndexBack[deep3dot][2] = true;
      }
    }
  }
  return;
}

inline void CopyNum(unsigned int hang, unsigned int* shu, char fuhao) {
  const auto& dic = dict[*shu];
  memcpy(CircleAnswerPtr[NowBuff][hang] + CircleAnswerIdx[NowBuff][hang],
         dic.charNum + dic.start, dic.length);
  CircleAnswerIdx[NowBuff][hang] += dic.length;
  CircleAnswerPtr[NowBuff][hang][CircleAnswerIdx[NowBuff][hang]++] = fuhao;
}

void Forward(register unsigned int NowRoot) {
  unsigned int deep1dot = NowRoot;  //第一个点
  IndexForw[deep1dot] = true;       //把当前点标记为路径上的点
  unsigned int ccount1 = ChildrenNum[deep1dot];
  for (unsigned int i = 0; i < ccount1; i++) {
    unsigned int Money1 = (ChildrenStart[deep1dot] + i)->second;
    unsigned int deep2dot = (ChildrenStart[deep1dot] + i)->first;  //第二个点
    if (IndexForw[deep2dot]) continue;
    IndexForw[deep2dot] = true;
    unsigned int ccount2 = ChildrenNum[deep2dot];
    for (unsigned int j = 0; j < ccount2; j++) {
      unsigned int Money2 = (ChildrenStart[deep2dot] + j)->second;
      if (Money2 * 5 < Money1 || Money2 > 3 * Money1) continue;
      unsigned int deep3dot = (ChildrenStart[deep2dot] + j)->first;  //第三个点
      if (IndexForw[deep3dot]) continue;
      IndexForw[deep3dot] = true;
      unsigned int ccount3 = ChildrenNum[deep3dot];
      for (unsigned int k = 0; k < ccount3; k++) {
        unsigned int Money3 = (ChildrenStart[deep3dot] + k)->second;
        if (Money3 * 5 < Money2 || Money3 > 3 * Money2) continue;
        unsigned int deep4dot =
            (ChildrenStart[deep3dot] + k)->first;  //第四个点
        if (deep4dot == NowRoot) {
          if ((Money1 * 5) < Money3 || Money1 > 3 * Money3) continue;
          CopyNum(0, &deep1dot, ',');
          CopyNum(0, &deep2dot, ',');
          CopyNum(0, &deep3dot, '\n');
          ++circleTimes;  //存三环
          continue;
        }
        if (IndexForw[deep4dot]) continue;
        IndexForw[deep4dot] = true;
        unsigned int ccount4 = ChildrenNum[deep4dot];
        for (unsigned int m = 0; m < ccount4; m++) {
          unsigned int Money4 = (ChildrenStart[deep4dot] + m)->second;
          if (Money4 * 5 < Money3 || Money4 > 3 * Money3) continue;
          unsigned int deep5dot =
              (ChildrenStart[deep4dot] + m)->first;  //第五个点
          if (deep5dot == NowRoot) {
            if ((Money1 * 5) < Money4 || Money1 > 3 * Money4) continue;
            CopyNum(1, &deep1dot, ',');
            CopyNum(1, &deep2dot, ',');
            CopyNum(1, &deep3dot, ',');
            CopyNum(1, &deep4dot, '\n');
            ++circleTimes;  //存四环
            continue;
          }
          if (!IndexBack[deep5dot][2] || IndexForw[deep5dot]) continue;
          IndexForw[deep5dot] = true;
          unsigned int ccount5 = ChildrenNum[deep5dot];
          for (unsigned int n = 0; n < ccount5; n++) {
            unsigned int Money5 = (ChildrenStart[deep5dot] + n)->second;
            if (Money5 * 5 < Money4 || Money5 > 3 * Money4) continue;
            unsigned int deep6dot =
                (ChildrenStart[deep5dot] + n)->first;  //第六个点
            if (deep6dot == NowRoot) {
              if ((Money1 * 5) < Money5 || Money1 > 3 * Money5) continue;
              CopyNum(2, &deep1dot, ',');
              CopyNum(2, &deep2dot, ',');
              CopyNum(2, &deep3dot, ',');
              CopyNum(2, &deep4dot, ',');
              CopyNum(2, &deep5dot, '\n');
              ++circleTimes;  //存五环
              continue;
            }
            if (!IndexBack[deep6dot][1] || IndexForw[deep6dot]) continue;
            IndexForw[deep6dot] = true;
            unsigned int ccount6 = ChildrenNum[deep6dot];
            for (unsigned int p = 0; p < ccount6; p++) {
              unsigned int Money6 = (ChildrenStart[deep6dot] + p)->second;
              if (Money6 * 5 < Money5 || Money6 > 3 * Money5) continue;
              unsigned int deep7dot =
                  (ChildrenStart[deep6dot] + p)->first;  //第七个点
              // unsigned int Money7 = (ChildrenStart[deep7dot]+p)->second;
              if (deep7dot == NowRoot) {
                if ((Money1 * 5) < Money6 || Money1 > 3 * Money6) continue;
                CopyNum(3, &deep1dot, ',');
                CopyNum(3, &deep2dot, ',');
                CopyNum(3, &deep3dot, ',');
                CopyNum(3, &deep4dot, ',');
                CopyNum(3, &deep5dot, ',');
                CopyNum(3, &deep6dot, '\n');
                ++circleTimes;  //存六环
                continue;
              }
              if (!IndexBack[deep7dot][0] || IndexForw[deep7dot]) {
                continue;
              } else {
                unsigned int ccount7 = ChildrenNum[deep7dot];
                for (unsigned int q = 0; q < ccount7; q++) {
                  unsigned int Money7 = (ChildrenStart[deep7dot] + q)->second;
                  if (Money7 * 5 < Money6 || Money7 > 3 * Money6) continue;
                  unsigned int deep8dot =
                      (ChildrenStart[deep7dot] + q)->first;  //第七个点
                  if (deep8dot == NowRoot) {
                    if ((Money1 * 5) < Money7 || Money1 > 3 * Money7) continue;
                    CopyNum(4, &deep1dot, ',');
                    CopyNum(4, &deep2dot, ',');
                    CopyNum(4, &deep3dot, ',');
                    CopyNum(4, &deep4dot, ',');
                    CopyNum(4, &deep5dot, ',');
                    CopyNum(4, &deep6dot, ',');
                    CopyNum(4, &deep7dot, '\n');
                    ++circleTimes;  //存七环
                    continue;
                  }
                }
              }
            }
            IndexForw[deep6dot] = false;
          }
          IndexForw[deep5dot] = false;
        }
        IndexForw[deep4dot] = false;
      }
      IndexForw[deep3dot] = false;
    }
    IndexForw[deep2dot] = false;
  }
  return;
}

void BothWayFind(unsigned int maxID, unsigned int ChildrenNum,
                 unsigned int* SyncPtr) {
  unsigned int LastEndID = 0;
  unsigned int NowStartID = 0;
  unsigned int NowEndID = 0;
  for (unsigned int i = 0; i < 5; i++) {
    CircleAnswerPtr[0][i] = CircleAnswer[i];
  }
  while (SyncPtr[5] != 1) {
    while (SyncPtr[4] == 1) {
      usleep(1);
    }
    SyncPtr[4] = 1;
    for (unsigned int i = 0; i < CHILDRENCORE + 1; i++) {
      if (SyncPtr[i] > NowStartID) NowStartID = SyncPtr[i];
    }
    NowEndID = NowStartID + maxID / 32;
    SyncPtr[ChildrenNum] = NowEndID;
    if (NowEndID >= maxID) {
      SyncPtr[5] = 1;
      NowEndID = maxID;
    }
    SyncPtr[4] = 0;
    CircleAnswerID[NowBuff][0] = NowStartID;
    CircleAnswerID[NowBuff][1] = NowEndID;
    // cout<<NowStartID<<","<<NowEndID<<endl;
    for (unsigned int i = LastEndID; i < NowStartID; i++)  //预处理
    {
      IndexForw[i] = true;
    }
    LastEndID = NowEndID;
    for (unsigned int i = NowStartID; i < NowEndID; i++)  //双向dfs
    {
      if (!IndexForw[i]) {
        Count = 0;
        Back(i);
        Forward(i);
        for (unsigned int j = 0; j < Count; j++) {
          IndexBack[Flag[j]][0] = false;
          IndexBack[Flag[j]][1] = false;
          IndexBack[Flag[j]][2] = false;
        }
      }
    }
    NowBuff++;
    if (NowBuff == 9) break;
    for (unsigned int i = 0; i < 5; i++) {
      CircleAnswerPtr[NowBuff][i] =
          &CircleAnswerPtr[NowBuff - 1][i][CircleAnswerIdx[NowBuff - 1][i]];
    }
  }
}

int main() {
#ifdef TIME
  cout << "开始加载邻接图.....<<endl";
  steady_clock::time_point startLoad = steady_clock::now();
#endif

  unsigned int maxID = LoadMap();

#ifdef TIME
  cout << "节点数：" << maxID << ",时间："
       << duration_cast<duration<double>>(steady_clock::now() - startLoad)
              .count()
       << ",开始启动子进程....." << endl;
  steady_clock::time_point startBuild = steady_clock::now();
#endif

  pid_t pid_Children[CHILDRENCORE];
  unsigned int ChildrenNum = 0;
  unsigned int CircleTimesShm = shmget(
      IPC_PRIVATE, (CHILDRENCORE + 1) * sizeof(unsigned int), IPC_CREAT | 0600);
  unsigned int OKShm = shmget(
      IPC_PRIVATE, (CHILDRENCORE + 1) * sizeof(unsigned int), IPC_CREAT | 0600);
  unsigned int OUTShm =
      shmget(IPC_PRIVATE, sizeof(unsigned int), IPC_CREAT | 0600);
  unsigned int SyncShm =
      shmget(IPC_PRIVATE, 6 * sizeof(unsigned int), IPC_CREAT | 0600);
  for (unsigned int i = 0; i < CHILDRENCORE; i++) {
    pid_Children[i] = 0;
  }
  for (unsigned int i = 0; i < CHILDRENCORE; i++) {
    if (ChildrenNum == 0) {
      pid_Children[i] = fork();
    }
    if (pid_Children[i] == 0) {
      ChildrenNum = i + 1;
      break;
    }
  }
  unsigned int* CircleTimesPtr = (unsigned int*)shmat(CircleTimesShm, NULL, 0);
  unsigned int* OKPtr = (unsigned int*)shmat(OKShm, NULL, 0);
  unsigned int* OUTPtr = (unsigned int*)shmat(OUTShm, NULL, 0);
  unsigned int* SyncPtr = (unsigned int*)shmat(SyncShm, NULL, 0);
  *OUTPtr = -1;

#ifdef TIME
  cout << "子进程" << ChildrenNum << "启动成功,时间："
       << duration_cast<duration<double>>(steady_clock::now() - startBuild)
              .count()
       << ",开始找环....." << endl;
  steady_clock::time_point startFind = steady_clock::now();
#endif

  BothWayFind(maxID, ChildrenNum, SyncPtr);
  CircleTimesPtr[ChildrenNum] = circleTimes;
  OKPtr[ChildrenNum] = -1;
  OKPtr[ChildrenNum] = -1;
  OKPtr[ChildrenNum] = -1;

#ifdef TIME
  cout << "进程：" << ChildrenNum << "找到" << circleTimes << "个环,时间：";
  cout << duration_cast<duration<double>>(steady_clock::now() - startFind)
              .count()
       << endl;
  steady_clock::time_point startWait = steady_clock::now();
#endif
  while (true) {
    unsigned int i = 0;
    for (i = 0; i < CHILDRENCORE + 1; i++) {
      if (OKPtr[i] != -1) break;
    }
    if (i == CHILDRENCORE + 1) {
      break;
    }
    usleep(1);
  }
  unsigned int AllCircle = 0;
  if (ChildrenNum == 0) {
    for (unsigned int i = 0; i < CHILDRENCORE + 1; i++) {
      AllCircle += CircleTimesPtr[i];
    }
  }
#ifdef TIME
  if (ChildrenNum == 0) cout << "总环数：" << AllCircle << endl;
  steady_clock::time_point startOut = steady_clock::now();
#endif

  SaveAnswer(maxID, AllCircle, ChildrenNum, OUTPtr);

#ifdef TIME
  cout << "输出结束,时间：";
  cout
      << duration_cast<duration<double>>(steady_clock::now() - startOut).count()
      << endl;
#endif
  if (ChildrenNum != 0) {
    exit(0);
  }
  return 0;
}
