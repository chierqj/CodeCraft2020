// #include <arm_neon.h>
#include <fcntl.h>
#include <sys/ipc.h>
#include <sys/mman.h>
#include <sys/shm.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#ifdef LOCAL
const std::string LOCAL_TRAIN = "../data/train_data.txt";
const std::string LOCAL_PREDICT = "../data/test_data.txt";
const std::string LOCAL_RESULT = "../data/result.txt";
const std::string LOCAL_ANSWER = "../data/answer.txt";
#else
const std::string TRAIN = "/data/train_data.txt";
const std::string PREDICT = "/data/test_data.txt";
const std::string RESULT = "/projects/student/result.txt";
const std::string ANSWER = "/projects/student/answer.txt";
#endif

const int NTHREAD = 4;
bool FLAG_TRAIN[NTHREAD] = {false};
bool FLAG_PREDICT[NTHREAD] = {false};

void Train(int pid) {
  // std::cerr << "pid: " << pid << "\n";
  sleep(1);
}
void Predict(int pid) {}

int main() {
  pid_t child1 = 0, child2 = 0, child3 = 0;
  child1 = fork();
  if (child1) child2 = fork();
  if (child2) child3 = fork();

  int num = 0;
  if (child1 == 0) {
    num = 1;
  } else if (child2 == 0) {
    num = 2;
  } else if (child3 == 0) {
    num = 3;
  } else {
    num = 0;
  }
  std::cerr << num << "\n";
  Train(num);
  while (!FLAG_TRAIN[0]) {
    sleep(1);
  }
  return 0;
}