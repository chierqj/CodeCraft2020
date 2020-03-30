#include <bits/stdc++.h>
#include <fcntl.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
using namespace std;

#define columnLen 1000
#define PROCESS_NUMBER 8
#define THREAD_NUMBER 8

// 线上参数
const int trainL = 1800;
const int hollowNight = 4000;
const int smallSet = hollowNight / PROCESS_NUMBER;

string trainFile = "/data/train_data.txt";
string testFile = "/data/test_data.txt";
string predictFile = "/projects/student/result.txt";

const int columnSize = columnLen * sizeof(int);
const float ZeroRate = 0.047;

int testFileSize, totalRows, rowsPerThread;
char *trainBuf, *testBuf, *outputBuf;

int param[1024], param0[THREAD_NUMBER * 1024], param1[THREAD_NUMBER * 1024];
int trainCnt0[THREAD_NUMBER], trainCnt1[THREAD_NUMBER];

void magicTest(int id) {
  vector<int> record0, record1;
  int rows = id != PROCESS_NUMBER - 1
                 ? rowsPerThread
                 : totalRows - rowsPerThread * (PROCESS_NUMBER - 1);
  int prows = rows;
  size_t bufPtr = 6000 * rowsPerThread * id;
  char* buf = testBuf;

  int cnt = 0;

  while (rows--) {
    int score = 0, nxtPtr = bufPtr + 6000;

    if (cnt >= smallSet && buf[bufPtr + 2] <= '9') {
      record0.push_back(100000000);
      record1.push_back(100000000);
      bufPtr = nxtPtr;
      continue;
    }

    cnt += 1;
    for (int i = 0; i < columnLen; i += 4) {
      int p0 = buf[bufPtr + 2] - '0';
      int p1 = buf[bufPtr + 8] - '0';
      int p2 = buf[bufPtr + 14] - '0';
      int p3 = buf[bufPtr + 20] - '0';
      score = score + (param[i] * p0) + (param[i + 1] * p1) +
              (param[i + 2] * p2) + (param[i + 3] * p3);
      bufPtr += 24;
    }
    record0.push_back(score);
    record1.push_back(score);
  }
  // cout<<cnt<<endl;
  sort(record0.begin(), record0.end());
  int limitPos = prows * ZeroRate;
  int limit = record0[limitPos];

  int outputFd = open(predictFile.c_str(), O_RDWR | O_CREAT, 0777);
  outputBuf = (char*)mmap(nullptr, totalRows << 1, PROT_READ | PROT_WRITE,
                          MAP_SHARED, outputFd, 0);

  int idx = rowsPerThread * id;
  for (int i = idx; i < idx + prows; i++) {
    if (record1[i - idx] < limit) {
      outputBuf[i << 1] = '0';
    } else {
      outputBuf[i << 1] = '1';
    }
  }
  // munmap((void*)outputBuf, totalRows << 1);
  exit(0);
}

void* trainThread(void* arg) {
  int id = *(int*)arg, rowsPerThread = trainL / THREAD_NUMBER,
      rows = rowsPerThread;
  size_t bufPtr = 6002 * rowsPerThread * id, memOffset = id << 10;
  char* buf = trainBuf;

  int local[1024];

  while (bufPtr != 0 && buf[bufPtr - 1] != '\n') bufPtr++;

  while (rows--) {
    int nxtPtr = bufPtr + 6000, flag = 0;
    for (size_t i = 0; i < columnLen; i += 4) {
      if (buf[bufPtr + 19] == '.') {
        local[i] = buf[bufPtr + 2] - '0';
        local[i + 1] = buf[bufPtr + 8] - '0';
        local[i + 2] = buf[bufPtr + 14] - '0';
        local[i + 3] = buf[bufPtr + 20] - '0';
        bufPtr += 24;
      } else {
        if (buf[bufPtr] == '-') {
          local[i] = '0' - buf[bufPtr + 3];
          bufPtr += 7;
        } else {
          local[i] = buf[bufPtr + 2] - '0';
          bufPtr += 6;
        }
        i -= 3;
      }
    }

    // if (flag) {
    //     while (buf[bufPtr + 1] != '\n') bufPtr++;
    //     continue;
    // }

    if (buf[bufPtr] == '0') {
      trainCnt0[id] += 1;
      for (size_t i = 0; i < columnLen; i += 4) {
        param0[memOffset + i] += local[i];
        param0[memOffset + i + 1] += local[i + 1];
        param0[memOffset + i + 2] += local[i + 2];
        param0[memOffset + i + 3] += local[i + 3];
      }
    } else {
      trainCnt1[id] += 1;
      for (size_t i = 0; i < columnLen; i += 4) {
        param1[memOffset + i] += local[i];
        param1[memOffset + i + 1] += local[i + 1];
        param1[memOffset + i + 2] += local[i + 2];
        param1[memOffset + i + 3] += local[i + 3];
      }
    }

    bufPtr += 2;
  }
  pthread_exit(nullptr);
}

// 共享: param0, param1, cnt0, cnt1, param, recs
pthread_t myThread[THREAD_NUMBER];
int sParam0[1024], sParam1[1024];

int main() {
  int trainFd = open(trainFile.c_str(), O_RDONLY);
  int testFd = open(testFile.c_str(), O_RDONLY);
  int testFileSize = lseek(testFd, 0, SEEK_END);
  totalRows = testFileSize / 6000;
  rowsPerThread = totalRows / PROCESS_NUMBER;

  trainBuf =
      (char*)mmap(nullptr, 7010 * trainL, PROT_READ, MAP_SHARED, trainFd, 0);
  testBuf =
      (char*)mmap(nullptr, testFileSize, PROT_READ, MAP_SHARED, testFd, 0);
  // -------- out put info ----------
  int outputFileSize = totalRows << 1;
  char* tmp = new char[outputFileSize];
  memset(tmp, '\n', outputFileSize);

  FILE* fp = fopen(predictFile.c_str(), "w");
  fwrite(tmp, 1, outputFileSize, fp);
  fclose(fp);
  // --------- init finished ---------

  // --------- trainning -------------
  for (int i = 0; i < THREAD_NUMBER; i++) {
    int* x = new int(i);
    pthread_create(&myThread[i], nullptr, trainThread, (void*)x);
  }
  int cnt0 = 0, cnt1 = 0;
  for (int i = 0; i < THREAD_NUMBER; i++) {
    pthread_join(myThread[i], nullptr);
    size_t offset = i << 10;
    cnt0 += trainCnt0[i];
    cnt1 += trainCnt1[i];
    for (int j = 0; j < columnLen; j++) {
      sParam0[j] += param0[offset + j];
      sParam1[j] += param1[offset + j];
    }
  }
  // cout<<cnt0<<" "<<cnt1<<endl;
  for (int i = 0; i < columnLen; i++) {
    param[i] = (int)(1000.0 * sParam1[i] / cnt1 - 1000.0 * sParam0[i] / cnt0);
  }
  // ---------- test -------------------
  for (int i = 0; i < PROCESS_NUMBER; i++) {
    if (i == PROCESS_NUMBER - 1) {
      magicTest(i);
      break;
    }
    pid_t pid = fork();
    if (pid == 0) {
      magicTest(i);
    }
  }
}