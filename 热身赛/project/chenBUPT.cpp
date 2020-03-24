#include <arm_neon.h>
#include <fcntl.h>
#include <sys/ipc.h>
#include <sys/mman.h>
#include <sys/shm.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
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
/***********************************************/
//#define TEST
#ifdef TEST
steady_clock::time_point start = steady_clock::now();
#endif
/***********************************************/
#define M 7500
#define N 1000
#define NT 1000
#define TestNum 20000
char answer[50000];
int trainDataSet[M][N];
int testData[TestNum][N];
int one = 0;
int zero = 0;
int mean0[N];
int mean1[N];
int16_t meanSum[N];
int16_t meanSub[N];
#define FATHER 0
#define CHILDREN1 1
#define CHILDREN2 2
#define CHILDREN3 3

#define TEST1

#ifdef TEST
#define trainFlie "./data/train_data.txt"
#define testFile "./data/test_data.txt"
#define resultFile "./projects/student/result.txt"
string predictFile = "./projects/student/result.txt";
string answerFile = "./projects/student/answer.txt";
string gt[] = {"", "\t\t\t", "\t\t\t\t\t\t", "\t\t\t\t\t\t\t\t\t"};
#else
#define trainFlie "/data/train_data.txt"
#define testFile "/data/test_data.txt"
#define resultFile "/projects/student/result.txt"
string predictFile = "/projects/student/result.txt";
string answerFile = "/projects/student/answer.txt";
#endif

bool loadAnswerData(string awFile, vector<int> &awVec) {
  ifstream infile(awFile.c_str());
  while (infile) {
    string line;
    int aw;
    getline(infile, line);
    if (line.size() > 0) {
      stringstream sin(line);
      sin >> aw;
      awVec.push_back(aw);
    }
  }
  infile.close();
  return true;
}
int loadTestData(int serialNumber) {
  int Max[] = {5000, 5000, 5000, 5000};
  long i = 0;
  int8_t *p;
  int8_t *data = NULL;
  int fd = open(testFile, O_RDONLY);
  long size = lseek(fd, 0, SEEK_END);
  data = (int8_t *)mmap(NULL, size, PROT_READ, MAP_PRIVATE, fd, 0);
  long n = 0, m = 0;
  int8_t *quarter[4] = {
      data, data + (int)(size * (Max[0]) / TestNum),
      data + (int)(size * (Max[0] + Max[1]) / TestNum),
      data + (int)(size * (Max[0] + Max[1] + Max[2]) / TestNum)};
  while (*(quarter[1]) != '\n') quarter[1]++;
  quarter[1]++;
  while (*(quarter[2]) != '\n') quarter[2]++;
  quarter[2]++;
  while (*(quarter[3]) != '\n') quarter[3]++;
  quarter[3]++;
  int8_t *dataCopy = data;
  data = quarter[serialNumber];
  int distance = 0;

  int8x16_t Twenty = vdupq_n_s8(20);
  int8x16_t Ten = vdupq_n_s8(10);
  int8x16_t ZeroChar = vdupq_n_s8('0');
  int8x16_t One8 = vdupq_n_s8(0);
  int8x16_t Two8 = vdupq_n_s8(0);

  int16x8_t One16_Low = vdupq_n_s16(0);
  int16x8_t Two16_Low = vdupq_n_s16(0);
  int16x8_t One16_High = vdupq_n_s16(0);
  int16x8_t Two16_High = vdupq_n_s16(0);

  int32x4_t Sum_Low_Low = vdupq_n_s32(0);
  int32x4_t Sum_Low_High = vdupq_n_s32(0);
  int32x4_t Sum_High_Low = vdupq_n_s32(0);
  int32x4_t Sum_High_High = vdupq_n_s32(0);

  int16_t *meanSumPtr = meanSum;
  int16_t *meanSubPtr = meanSub;
  int Flag = 0;
  while (true) {
    p = data + i;
    if (n >= 62) {
      distance =
          vgetq_lane_s32(Sum_Low_Low, 0) + vgetq_lane_s32(Sum_Low_Low, 1) +
          vgetq_lane_s32(Sum_Low_Low, 2) + vgetq_lane_s32(Sum_Low_Low, 3) +
          vgetq_lane_s32(Sum_Low_High, 0) + vgetq_lane_s32(Sum_Low_High, 1) +
          vgetq_lane_s32(Sum_Low_High, 2) + vgetq_lane_s32(Sum_Low_High, 3) +
          vgetq_lane_s32(Sum_High_Low, 0) + vgetq_lane_s32(Sum_High_Low, 1) +
          vgetq_lane_s32(Sum_High_Low, 2) + vgetq_lane_s32(Sum_High_Low, 3) +
          vgetq_lane_s32(Sum_High_High, 0) + vgetq_lane_s32(Sum_High_High, 1) +
          vgetq_lane_s32(Sum_High_High, 2) + vgetq_lane_s32(Sum_High_High, 3);
      meanSumPtr = meanSum;
      meanSubPtr = meanSub;
      Sum_Low_Low = vdupq_n_s32(0);
      Sum_Low_High = vdupq_n_s32(0);
      Sum_High_Low = vdupq_n_s32(0);
      Sum_High_High = vdupq_n_s32(0);
      i += 48;
      p = data + i;
      answer[m] = distance >= 0 ? '0' : '1';
      distance = 0;
      m++;
      n = 0;
      if ((serialNumber < 3 && p >= quarter[serialNumber + 1]) ||
          (serialNumber == 3 && p - dataCopy >= size)) {
        return m;
        break;
      }
    }
    One8 = vld1q_lane_s8((p + 2), One8, 0);
    Two8 = vld1q_lane_s8((p + 3), Two8, 0);
    One8 = vld1q_lane_s8((p + 8), One8, 1);
    Two8 = vld1q_lane_s8((p + 9), Two8, 1);
    One8 = vld1q_lane_s8((p + 14), One8, 2);
    Two8 = vld1q_lane_s8((p + 15), Two8, 2);
    One8 = vld1q_lane_s8((p + 20), One8, 3);
    Two8 = vld1q_lane_s8((p + 21), Two8, 3);
    One8 = vld1q_lane_s8((p + 26), One8, 4);
    Two8 = vld1q_lane_s8((p + 27), Two8, 4);
    One8 = vld1q_lane_s8((p + 32), One8, 5);
    Two8 = vld1q_lane_s8((p + 33), Two8, 5);
    One8 = vld1q_lane_s8((p + 38), One8, 6);
    Two8 = vld1q_lane_s8((p + 39), Two8, 6);
    One8 = vld1q_lane_s8((p + 44), One8, 7);
    Two8 = vld1q_lane_s8((p + 45), Two8, 7);
    One8 = vld1q_lane_s8((p + 50), One8, 8);
    Two8 = vld1q_lane_s8((p + 51), Two8, 8);
    One8 = vld1q_lane_s8((p + 56), One8, 9);
    Two8 = vld1q_lane_s8((p + 57), Two8, 9);
    One8 = vld1q_lane_s8((p + 62), One8, 10);
    Two8 = vld1q_lane_s8((p + 63), Two8, 10);
    One8 = vld1q_lane_s8((p + 68), One8, 11);
    Two8 = vld1q_lane_s8((p + 69), Two8, 11);
    One8 = vld1q_lane_s8((p + 74), One8, 12);
    Two8 = vld1q_lane_s8((p + 75), Two8, 12);
    One8 = vld1q_lane_s8((p + 80), One8, 13);
    Two8 = vld1q_lane_s8((p + 81), Two8, 13);
    One8 = vld1q_lane_s8((p + 86), One8, 14);
    Two8 = vld1q_lane_s8((p + 87), Two8, 14);
    One8 = vld1q_lane_s8((p + 92), One8, 15);
    Two8 = vld1q_lane_s8((p + 93), Two8, 15);

    One8 = vsubq_s8(One8, ZeroChar);  //*(p+2) - '0'
    Two8 = vsubq_s8(Two8, ZeroChar);  //*(p+3) - '0'
    One8 = vmulq_s8(One8, Ten);       //(*(p+2) - '0')*10

    One16_Low = vmull_s8(vget_low_s8(One8), vget_low_s8(Twenty));
    One16_High = vmull_s8(vget_high_s8(One8),
                          vget_low_s8(Twenty));  // 200*(*(p+2) - '0')
    Two16_Low = vmull_s8(vget_low_s8(Two8), vget_low_s8(Twenty));
    Two16_High =
        vmull_s8(vget_high_s8(Two8), vget_low_s8(Twenty));  // 20*(*(p+3) - '0')

    One16_Low = vaddq_s16(One16_Low, Two16_Low);
    One16_High = vaddq_s16(
        One16_High, Two16_High);  // 200*(*(p+2) - '0') + 20*(*(p+3) - '0')

    Two16_Low = vld1q_s16(meanSumPtr);
    meanSumPtr += 8;
    Two16_High = vld1q_s16(meanSumPtr);
    meanSumPtr += 8;  // meanSum
    One16_Low = vsubq_s16(One16_Low, Two16_Low);
    One16_High = vsubq_s16(
        One16_High,
        Two16_High);  // 200*(*(p+2) - '0') + 20*(*(p+3) - '0') - meanSumPtr

    Two16_Low = vld1q_s16(meanSubPtr);
    meanSubPtr += 8;
    Two16_High = vld1q_s16(meanSubPtr);
    meanSubPtr += 8;  // meanSum
    One16_Low = vmulq_s16(One16_Low, Two16_Low);
    One16_High =
        vmulq_s16(One16_High,
                  Two16_High);  //(200*(*(p+2) - '0')... - meanSumPtr) * meanSub

    Sum_Low_Low = vaddw_s16(Sum_Low_Low, vget_low_s16(One16_Low));
    Sum_Low_High = vaddw_s16(Sum_Low_High, vget_high_s16(One16_Low));
    Sum_High_Low = vaddw_s16(Sum_High_Low, vget_low_s16(One16_High));
    Sum_High_High = vaddw_s16(Sum_High_High, vget_high_s16(One16_High));
    i += 96;
    n++;
  }
  return true;
}
void loadTrainData(int serialNumber, int *mean0Ptr, int *mean1Ptr, int *zeroPtr,
                   int *onePtr) {
  int Max[] = {M / 4, M / 4, M / 4, M / 4};
  int Now = 0;
  long n = N - 1, m = 0, i = 0;
  char *data = NULL, *p = NULL;
  int fd = open(trainFlie, O_RDONLY);
  long size = lseek(fd, 0, SEEK_END);
  data = (char *)mmap(NULL, size, PROT_READ, MAP_PRIVATE, fd, 0);
  data = data + size - 1;
  char *quarter[4] = {data, data - (int)(6002 * (Max[0])),
                      data - (int)(6002 * (Max[0] + Max[1])),
                      data - (int)(6002 * (Max[0] + Max[1] + Max[2]))};
  while (*(quarter[1]) != '\n') quarter[1]--;
  while (*(quarter[2]) != '\n') quarter[2]--;
  while (*(quarter[3]) != '\n') quarter[3]--;
  data = quarter[serialNumber];
  data--;
  Now = *(data) - '0';
  if (Now)
    one++;
  else
    zero++;
  i += 2;
#ifdef TEST
  cout << gt[serialNumber] << serialNumber << " Deal OK: ";
  cout << duration_cast<duration<double>>(steady_clock::now() - start).count()
       << endl;
#endif
  while (true) {
    p = data - i;
    if (n < 0) {
      m++;
      if ((serialNumber < 3 && p <= quarter[serialNumber + 1]) ||
          (serialNumber == 3 && m >= Max[serialNumber])) {
        break;
      } else {
        Now = *(p) - '0';
        if (Now)
          one++;
        else
          zero++;
      }
      n = N - 1;
      i += 2;
      p = data - i;
    }
    if (*(p - 5) == '-') {
      if (Now)
        mean1[n] += ('0' - *(p - 1)) * 10 + 100 * ('0' - *(p - 2));
      else
        mean0[n] += ('0' - *(p - 1)) * 10 + 100 * ('0' - *(p - 2));
      i += 7;
    } else {
      if (Now)
        mean1[n] += (*(p - 1) - '0') * 10 + 100 * (*(p - 2) - '0');
      else
        mean0[n] += (*(p - 1) - '0') * 10 + 100 * (*(p - 2) - '0');
      i += 6;
    }
    n--;
  }
  for (int i = 0; i < N; i++) {
    mean0Ptr[i] += mean0[i];
    mean1Ptr[i] += mean1[i];
  }
  *zeroPtr += zero;
  *onePtr += one;
}

int main() {
#ifdef TEST
  vector<int> answerVec;
  vector<int> predictVec;
  int correctCount;
  double accurate;
#endif

  /****************¹²ÏíµÄÄÚ´æ±äÁ¿*******************/
  pid_t pid_Children1, pid_Children2, pid_Children3;
  int resultShm[4];
  char *resultPtr[4];
  int mean0Shm = shmget(IPC_PRIVATE, N * sizeof(int), IPC_CREAT | 0600);
  int mean1Shm = shmget(IPC_PRIVATE, N * sizeof(int), IPC_CREAT | 0600);
  int zeroShm = shmget(IPC_PRIVATE, sizeof(int), IPC_CREAT | 0600);
  int oneShm = shmget(IPC_PRIVATE, sizeof(int), IPC_CREAT | 0600);
  int answerShm =
      shmget(IPC_PRIVATE, 2 * TestNum * sizeof(char), IPC_CREAT | 0600);
  int trainFlagShm = shmget(IPC_PRIVATE, 5 * sizeof(int), IPC_CREAT | 0600);
  int predictFlagShm = shmget(IPC_PRIVATE, 4 * sizeof(int), IPC_CREAT | 0600);
/****************¿ªÆô×Ó½ø³Ì£¨1-2-3£©*******************/
#ifdef TEST
  cout << "Farther" << gt[1] << "Children1"
       << "\t\t"
       << "Children2"
       << "\t\t"
       << "Children3" << endl;
#endif
  pid_Children1 = fork();
  if (pid_Children1) pid_Children2 = fork();
  if (pid_Children2) pid_Children3 = fork();
  int Num = 0;
  if (!pid_Children1)
    Num = CHILDREN1;
  else if (!pid_Children2)
    Num = CHILDREN2;
  else if (!pid_Children3)
    Num = CHILDREN3;
  else
    Num = FATHER;

#ifdef TEST
  cout << gt[Num] << Num << "  Start: ";
  cout << duration_cast<duration<double>>(steady_clock::now() - start).count()
       << endl;
#endif
  /****************¸÷×Ô»ñµÃ¾ä±ú*******************/
  int *mean0Ptr = (int *)shmat(mean0Shm, NULL, 0);
  int *mean1Ptr = (int *)shmat(mean1Shm, NULL, 0);
  int *zeroPtr = (int *)shmat(zeroShm, NULL, 0);
  int *onePtr = (int *)shmat(oneShm, NULL, 0);
  char *answerPtr = (char *)shmat(answerShm, NULL, 0);
  int *trainFlagPtr = (int *)shmat(trainFlagShm, NULL, 0);
  int *predictFlagPtr = (int *)shmat(predictFlagShm, NULL, 0);
  *trainFlagPtr = 0;
  /********************½ø³Ì¿ªÊ¼**************************/
  if (Num <= 3) {
/********¼ÓÔØÊý¾Ý**********/
#ifdef TEST
    cout << gt[Num] << Num << "  Loading Start: ";
    cout << duration_cast<duration<double>>(steady_clock::now() - start).count()
         << endl;
#endif
    loadTrainData(Num, mean0Ptr, mean1Ptr, zeroPtr, onePtr);
#ifdef TEST
    cout << gt[Num] << Num << "  Loading OK: ";
    cout << duration_cast<duration<double>>(steady_clock::now() - start).count()
         << endl;
#endif
    if (*zeroPtr + *onePtr > M - 100) {
      trainFlagPtr[4] = 1;
    }
    while ((trainFlagPtr[4]) == 0) {
      sleep(0.00001);
    }

    zero = *zeroPtr;
    one = *onePtr;
    for (int i = Num * 250; i < (Num + 1) * 250; i++) {
      mean0Ptr[i] /= zero;
      mean1Ptr[i] /= one;
    }
    trainFlagPtr[Num] = 1;
    trainFlagPtr[Num] = 1;
    while (trainFlagPtr[0] == 0 || trainFlagPtr[1] == 0 ||
           trainFlagPtr[2] == 0 || trainFlagPtr[3] == 0) {
      sleep(0.00001);
    }
    for (int i = 0; i < N; i++) {
      meanSub[i] = mean0Ptr[i] - mean1Ptr[i];
      meanSum[i] = mean0Ptr[i] + mean1Ptr[i];
    }
#ifdef TEST
    cout << gt[Num] << Num << " Predict Start ";
    cout << duration_cast<duration<double>>(steady_clock::now() - start).count()
         << endl;
#endif
    predictFlagPtr[Num] = loadTestData(Num);
#ifdef TEST
    cout << gt[Num] << Num << " Predict OK ";
    cout << duration_cast<duration<double>>(steady_clock::now() - start).count()
         << endl;
#endif
    while (!predictFlagPtr[0] || !predictFlagPtr[1] || !predictFlagPtr[2] ||
           !predictFlagPtr[3]) {
      sleep(0.00001);
    }
    int Start = 0;
    for (int i = 0; i < Num; i++) {
      Start += predictFlagPtr[i];
    }
    Start *= 2;
    for (int i = 0; i < predictFlagPtr[Num]; i++) {
      answerPtr[2 * i + Start] = answer[i];
      answerPtr[2 * i + Start + 1] = '\n';
    }
#ifdef TEST
    cout << gt[Num] << Num << " Exit ";
    cout << duration_cast<duration<double>>(steady_clock::now() - start).count()
         << endl;
#endif
    if (Num) exit(0);
  }
  waitpid(pid_Children1, NULL, 0);
  waitpid(pid_Children2, NULL, 0);
  waitpid(pid_Children3, NULL, 0);
  FILE *fp = fopen(resultFile, "wt");
  fputs(answerPtr, fp);
  fclose(fp);

#ifdef TEST
  cout << "All Time: ";
  cout << duration_cast<duration<double>>(steady_clock::now() - start).count()
       << endl;
  loadAnswerData(answerFile, answerVec);
  loadAnswerData(predictFile, predictVec);
  cout << "test data set size is " << predictVec.size() << endl;
  correctCount = 0;
  for (int j = 0; j < predictVec.size(); j++) {
    if (j < answerVec.size()) {
      if (answerVec[j] == predictVec[j]) {
        correctCount++;
      }
    } else {
      cout << "answer size less than the real predicted value" << endl;
    }
  }
  accurate = ((double)correctCount) / answerVec.size();
  cout << "the prediction accuracy is " << accurate << endl;
#endif
  return 0;
}
