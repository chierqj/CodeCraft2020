/**
 *    author:  ddd
 *    created: 2020.03.11
 **/

#include <bits/stdc++.h>
#include <fcntl.h>
#include <pthread.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#define unlikely(x) __builtin_expect(!!(x), 0)
#define is_digit(c) (c >= '0')

using namespace std;

// 读入训练集的量。
const int BUFSIZE_TRAIN = 11 << 20;

const int THREAD_COUNT = 8;
const int PROC_COUNT = 8;
char* buf;

vector<int> cof;

int test_data_size;

struct TrainInfo {
  int l;
  int r;
  vector<int> ans;
};

void* magic_train(void* train_info) {
  TrainInfo* info = (TrainInfo*)train_info;
  int l = info->l, r = info->r;
  vector<int> cur_cof(1000);
  vector<int> cnt(2);
  int j, ptr = 0;
  char* magic_buf = buf + l;
  const char* buf_end = buf + r;
  while (magic_buf < buf_end) {
    if (*(magic_buf + 6001) == '\n') {
#define go(op)                            \
  for (int i = 0; i < 250; ++i) {         \
    cur_cof[ptr] op*(magic_buf + 2);      \
    cur_cof[ptr + 1] op*(magic_buf + 8);  \
    cur_cof[ptr + 2] op*(magic_buf + 14); \
    cur_cof[ptr + 3] op*(magic_buf + 20); \
    magic_buf += 24;                      \
    ptr += 4;                             \
  }
      if (*(magic_buf + 6000) == '1') {
        go(+=);
        ++cnt[1];
      } else {
        go(-=);
        ++cnt[0];
      }
      magic_buf += 2;
      ptr = 0;
    } else {
      magic_buf = (char*)memchr(magic_buf + 6001, '\n', 1000) + 1;
    }
  }
  cur_cof.push_back(cnt[0]);
  cur_cof.push_back(cnt[0] + cnt[1]);
  info->ans = cur_cof;
  pthread_exit(NULL);
}

void train(const string& file_name) {
  int fd = open(file_name.c_str(), O_RDONLY);
  int sz = BUFSIZE_TRAIN;
  buf = (char*)mmap(NULL, sz, PROT_READ, MAP_PRIVATE, fd, 0);
  vector<int> cur(1000);
  vector<int> cur_cof(1002);
  vector<int> cnt(2);
  const int S = sz;
  pthread_attr_t attr;
  void* status;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  pthread_t tids[THREAD_COUNT];
  TrainInfo infos[THREAD_COUNT];
  int last_p = 0;
  for (int tc = 0; tc < THREAD_COUNT; ++tc) {
    int p = (tc + 1) * S / THREAD_COUNT;
    while (buf[p - 1] != '\n') --p;
    infos[tc].l = last_p;
    infos[tc].r = p;
    last_p = p;
    int rc = pthread_create(&tids[tc], NULL, magic_train, (void*)&(infos[tc]));
    if (rc) {
      cerr << "Wrong!" << endl;
      exit(-1);
    }
  }
  pthread_attr_destroy(&attr);
  for (int tc = 0; tc < THREAD_COUNT; ++tc) {
    int rc = pthread_join(tids[tc], &status);
    if (rc) {
      cerr << "Wrong!" << endl;
      exit(-1);
    }
  }
  for (int tc = 0; tc < THREAD_COUNT; ++tc) {
    for (int i = 0; i < (int)infos[tc].ans.size(); i++) {
      cur_cof[i] += infos[tc].ans[i];
    }
  }
  cof = cur_cof;
}

int max_cof = INT_MIN;
int max_cof_p = -1;
int big_cnt[PROC_COUNT];
vector<int> magic_vec[PROC_COUNT];

void* magic_fun(void* thread_id) {
  int tid = *((int*)thread_id);
  int sum = 0;
  int ptr = 0;
  const int L = test_data_size * tid / PROC_COUNT;
  const int R = L + test_data_size / PROC_COUNT;
  char* magic_buf = buf + L;
  const char* buf_end = buf + R;
  int cur_cof[cof.size()];
  for (int i = 0; i < (int)cof.size(); ++i) cur_cof[i] = cof[i];
  int big = 0;
  while (magic_buf < buf_end) {
    if (*(magic_buf + 2 + max_cof_p * 6) * 10 +
            *(magic_buf + 3 + max_cof_p * 6) >=
        49 * 10 + 48 + 5) {
      sum = INT_MAX;
      magic_buf += 6000;
      ++big;
    } else {
      int v[8] = {0};
      for (int i = 0; i < 125; ++i) {
        __builtin_prefetch(cur_cof + ptr + 7, 0);
        __builtin_prefetch(magic_buf + 44, 0);
        v[0] += cur_cof[ptr] * (*(magic_buf + 2));
        v[1] += cur_cof[ptr + 1] * (*(magic_buf + 8));
        v[2] += cur_cof[ptr + 2] * (*(magic_buf + 14));
        v[3] += cur_cof[ptr + 3] * (*(magic_buf + 20));
        v[4] += cur_cof[ptr + 4] * (*(magic_buf + 26));
        v[5] += cur_cof[ptr + 5] * (*(magic_buf + 32));
        v[6] += cur_cof[ptr + 6] * (*(magic_buf + 38));
        v[7] += cur_cof[ptr + 7] * (*(magic_buf + 44));
        magic_buf += 48;
        ptr += 8;
      }
      for (int i = 0; i < 8; ++i) sum += v[i];
    }
    magic_vec[tid].push_back(sum);
    sum = 0;
    ptr = 0;
  }
  big_cnt[tid] = big;
}

void work(const string& file_name, const string& result_file, int pid) {
  int fd = open(file_name.c_str(), O_RDONLY);
  test_data_size = lseek(fd, 0, SEEK_END);
  buf = (char*)mmap(NULL, test_data_size, PROT_READ, MAP_PRIVATE, fd, 0);
  magic_fun((void*)&pid);
  vector<int> vec = magic_vec[pid];
  int v1 = cof[(int)cof.size() - 2];
  int v2 = cof[(int)cof.size() - 1];
  int tot = (int)vec.size() - big_cnt[pid];
  int p = tot * 0.32;
  vector<int> vec2 = vec;
  nth_element(vec2.begin(), vec2.begin() + p, vec2.end());
  int magic = vec2[p];
  int fd2 = open(result_file.c_str(), O_RDWR | O_CREAT, 00777);
  char* ans = (char*)mmap(NULL, test_data_size / 3000, PROT_READ | PROT_WRITE,
                          MAP_SHARED, fd2, 0);
  close(fd2);
  const int L = pid * test_data_size / 6000 / PROC_COUNT;
  const int R = L + (int)vec.size();
  for (int i = L; i < R; ++i) {
    ans[i << 1] = (vec[i - L] <= magic ? '0' : '1');
    ans[i << 1 | 1] = '\n';
  }
  munmap(ans, test_data_size / 3000);
}

double check_ans(const string& answer, const string& result) {
  FILE* f1 = fopen(answer.c_str(), "r");
  FILE* f2 = fopen(result.c_str(), "r");
  assert(f1 && f2);
  int correct = 0;
  int total = 0;
  int v1, v2;
  vector<int> cnt(2);
  while (fscanf(f1, "%d", &v1) != -1) {
    assert(fscanf(f2, "%d", &v2) != -1);
    total += 1;
    ++cnt[v2];
    if (v1 == v2) {
      correct += 1;
    } else {
      // printf("%d %d\n", v1, v2);
    }
  }
  fclose(f1);
  fclose(f2);
  return 1.0 * correct / total;
}

int gcd(int x, int y) { return y ? gcd(y, x % y) : x; }

void normalize() {
  int sum_cof = 0;
  for (int i = 0; i < (int)cof.size() - 2; ++i) {
    sum_cof += cof[i];
    if (cof[i] > max_cof) {
      max_cof = cof[i];
      max_cof_p = i;
    }
  }
  int g = 0;
  for (int i = 0; i < (int)cof.size() - 2; ++i) {
    cof[i] = cof[i] * ((int)cof.size() - 2) - sum_cof;
    g = gcd(abs(cof[i]), g);
  }
  for (int i = 0; i < (int)cof.size() - 2; ++i) {
    cof[i] /= g;
  }
}

int main() {
  string trainFile = "/data/train_data.txt";
  string testFile = "/data/test_data.txt";
  string predictFile = "/projects/student/result.txt";
  string answerFile = "";

#ifdef TESTSPEED
  trainFile = "./train_data.txt";
  predictFile = "./result.txt";
  testFile = "./ddd.txt";
  answerFile = "./ddd_answer.txt";
  auto start = std::chrono::system_clock::now();
#endif

  train(trainFile);
  normalize();

  const int ANS_SIZE = 20000;

  int fd = open(predictFile.c_str(), O_RDWR | O_CREAT, 00777);
  lseek(fd, ANS_SIZE * 2 - 1, SEEK_SET);
  write(fd, " ", 1);
  close(fd);

  int id;
  pid_t pid;
  vector<pid_t> pids;
  for (int i = 1; i < PROC_COUNT; i++) {
    id = i;
    pid = fork();
    if (pid <= 0) break;
    pids.push_back(pid);
  }

  if (pid == -1) {
    cerr << "bad" << endl;
  } else {
    if (pid == 0) {
      work(testFile, predictFile, id);
      exit(0);
    } else {
      work(testFile, predictFile, 0);
    }
  }

#ifdef TESTSPEED
  std::chrono::duration<double, std::milli> duration =
      std::chrono::system_clock::now() - start;
  std::cout << "time elapsed: " << duration.count() << "ms" << std::endl;
  printf("%.9f\n", check_ans(answerFile, predictFile));
#endif

  exit(0);
}