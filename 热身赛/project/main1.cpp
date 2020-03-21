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
/***********************************************/
//#define TEST
#define trainFlie "/data/train_data.txt"
#define testFile "/data/test_data.txt"
#define resultFile "/projects/student/result.txt"
/***********************************************/
#define M 7200
#define N 1000
#define NT 1000
#define TestNum 20000
char answer[50000];
int disposeCount = 0;
int trainDataSet[M][N];
int trainDataSet_lab[M];
int testData[TestNum][N];

#ifdef TEST
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
bool loadTestData(int batch, char *data, long size)  // 0.6s-2000
{
  long i = 0;
  char *p;
  long n = 0, m = 0;
  long third = size / 3;
  long twoThird = size * 2 / 3;
  char *dataCopy = data;
  while (*(data + third) != '\n') {
    third++;
  }
  third++;
  while (*(data + twoThird) != '\n') {
    twoThird++;
  }
  twoThird++;
  if (batch == 1) {
    data = data + third;
  }
  if (batch == 2) {
    data = data + twoThird;
  }
  while (true) {
    p = data + i;
    if (n >= N) {
      m++;
      disposeCount++;
      n = 0;
      if ((batch == 0 && (p - dataCopy) >= third) ||
          (batch == 1 && (p - dataCopy) >= twoThird) ||
          (batch == 2 && ((p)-dataCopy) >= size)) {
        break;
      }
    }

    if (*p == '-') {
      testData[m][n] = ('0' - *(p + 3)) * 100 + 10 * ('0' - *(p + 4));
      i += 7;
    } else {
      testData[m][n] = (*(p + 2) - '0') * 100 + 10 * (*(p + 3) - '0');
      i += 6;
    }
    n++;
  }
  return true;
}
int one = 0;
int zero = 0;
void loadTrainData(int *mean0Ptr, int *mean1Ptr) {
  int Now = 0;
  long n = 0, m = 0;
  long i = 0;
  char *p;
  char *data = NULL;
  int fd = open(trainFlie, O_RDONLY);
  long size = lseek(fd, 0, SEEK_END);
  data = (char *)mmap(NULL, size, PROT_READ, MAP_PRIVATE, fd, 0);
  data = data + size - 2;
  Now = *(data) - '0';
  if (Now) {
    one++;
  } else {
    zero++;
  }
  i += 2;
  while (i < size) {
    p = data - i;
    if (n >= N) {
      m++;
      Now = *(p) - '0';
      if (Now) {
        one++;
      } else {
        zero++;
      }
      n = 0;
      i += 2;
      p = data - i;
      if (m >= M) {
        break;
      }
    }

    if (*(p - 5) == '-') {
      if (Now)
        mean1Ptr[N - 1 - n] += ('0' - *(p - 1)) * 10 + 100 * ('0' - *(p - 2));
      else
        mean0Ptr[N - 1 - n] += ('0' - *(p - 1)) * 10 + 100 * ('0' - *(p - 2));
      i += 7;
    } else {
      if (Now)
        mean1Ptr[N - 1 - n] += (*(p - 1) - '0') * 10 + 100 * (*(p - 2) - '0');
      else
        mean0Ptr[N - 1 - n] += (*(p - 1) - '0') * 10 + 100 * (*(p - 2) - '0');
      i += 6;
    }
    n++;
  }
  for (int i = 0; i < N; i++) {
    mean1Ptr[i] /= one;
    mean0Ptr[i] /= zero;
  }
}

int main() {
#ifdef TEST
  time_t clk[10];
  double times;
  clk[0] = clock();
  vector<int> answerVec;
  vector<int> predictVec;
  int correctCount;
  double accurate;
#endif

  pid_t pid_loadTrain, pid_loadTest1, pid_loadTest2;
  int mean0Shm = shmget(IPC_PRIVATE, N * sizeof(int), IPC_CREAT | 0600);
  int mean1Shm = shmget(IPC_PRIVATE, N * sizeof(int), IPC_CREAT | 0600);

  pid_loadTrain = fork();

  int *mean0Ptr = (int *)shmat(mean0Shm, NULL, 0);
  int *mean1Ptr = (int *)shmat(mean1Shm, NULL, 0);
  int i = 0, j = 0, k = 0, l = 0;

  if (!pid_loadTrain) {
    /***********************children0*********************************/
    /*loadTrainData-start*/
    loadTrainData(mean0Ptr, mean1Ptr);
/*loadTrainData-end*/
#ifdef TEST
    clk[1] = clock();
    cout << "children0 " << double(clk[1] - clk[0]) / CLOCKS_PER_SEC << endl;
#endif
    exit(0);
  } else {
    char *data = NULL;
    int fd = open(testFile, O_RDONLY);
    long size = lseek(fd, 0, SEEK_END);
    data = (char *)mmap(NULL, size, PROT_READ, MAP_PRIVATE, fd, 0);

    int Children1[2];
    int Children2[2];
    pipe(Children1);
    pipe(Children2);
    int result1Shm =
        shmget(IPC_PRIVATE, 2 * TestNum * sizeof(char), IPC_CREAT | 0600);
    int FlagShm = shmget(IPC_PRIVATE, sizeof(int), IPC_CREAT | 0600);

    pid_loadTest1 = fork();

    char *result1Ptr = (char *)shmat(result1Shm, NULL, 0);
    int *FlagPtr = (int *)shmat(FlagShm, NULL, 0);
    if (!pid_loadTest1) {
      /***********************children1*********************************/
      loadTestData(1, data, size);
      int Delay = 0;
      close(Children1[1]);
      char line[3] = "NO";
      while (true) {
        int n = read(Children1[0], line, 3);
        if (n == 2) break;
      }
      int distance = 0;
      int j;
      for (j = 0; j < disposeCount; j++) {
        distance = 0;

        for (int i = 0; i < N; i++) {
          distance +=
              (testData[j][i] + testData[j][i] - mean0Ptr[i] - mean1Ptr[i]) *
              (mean0Ptr[i] - mean1Ptr[i]);
        }
        result1Ptr[2 * j] = distance >= 0 ? '0' : '1';
        result1Ptr[2 * j + 1] = '\n';
      }
      result1Ptr[2 * disposeCount] = 0;
#ifdef TEST
      clk[1] = clock();
      cout << "children1 " << double(clk[1] - clk[0]) / CLOCKS_PER_SEC << endl;
#endif
      exit(0);
    } else {
      int result2Shm =
          shmget(IPC_PRIVATE, 2 * TestNum * sizeof(char), IPC_CREAT | 0600);
      /***********************children2*********************************/

      pid_loadTest2 = fork();

      char *result2Ptr = (char *)shmat(result2Shm, NULL, 0);
      if (!pid_loadTest2) {
        loadTestData(2, data, size);
        int Delay = 0;
        close(Children2[1]);
        char line[2];
        while (true) {
          int n = read(Children2[0], line, 2);
          if (n == 2) break;
        }

        int distance = 0;
        int j;
        for (j = 0; j < disposeCount; j++) {
          distance = 0;
          for (int i = 0; i < N; i++) {
            distance +=
                (testData[j][i] + testData[j][i] - mean0Ptr[i] - mean1Ptr[i]) *
                (mean0Ptr[i] - mean1Ptr[i]);
          }
          result2Ptr[2 * j] = distance >= 0 ? '0' : '1';
          result2Ptr[2 * j + 1] = '\n';
        }
        result2Ptr[2 * disposeCount] = 0;
#ifdef TEST
        clk[1] = clock();
        cout << "children2 " << double(clk[1] - clk[0]) / CLOCKS_PER_SEC
             << endl;
#endif
        exit(0);
      } else {
        loadTestData(0, data, size);

        waitpid(pid_loadTrain, NULL, 0);
        close(Children1[0]);
        close(Children2[0]);
        write(Children2[1], "O", 2);
        write(Children1[1], "O", 2);

        int distance = 0;
        FILE *fp = fopen(resultFile, "wt");
        int j = 0;
        for (j = 0; j < disposeCount; j++) {
          distance = 0;
          for (int i = 0; i < N; i++) {
            distance +=
                (testData[j][i] + testData[j][i] - mean0Ptr[i] - mean1Ptr[i]) *
                (mean0Ptr[i] - mean1Ptr[i]);
          }
          answer[2 * j] = distance >= 0 ? '0' : '1';
          answer[2 * j + 1] = '\n';
        }
        waitpid(pid_loadTest1, NULL, 0);
        j = 0;
        while (true) {
          if (result1Ptr[j] == 0) break;
          answer[2 * disposeCount + j] = result1Ptr[j];
          j++;
        }
        int k = 0;
        waitpid(pid_loadTest2, NULL, 0);
        while (true) {
          answer[2 * disposeCount + j] = result2Ptr[k];
          if (result2Ptr[k] == 0) break;
          k++;
          j++;
        }
        fputs(answer, fp);
        fclose(fp);
#ifdef TEST
        clk[1] = clock();
        cout << "parents" << double(clk[1] - clk[0]) / CLOCKS_PER_SEC << endl;
#endif
      }
    }
#ifdef TEST
    clk[5] = clock();
#endif

#ifdef TEST
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
  }
#ifdef TEST
  times = double(clk[5] - clk[0]) / CLOCKS_PER_SEC;
  cout << "all:" << times << endl;
#endif
  return 0;
}
