#include <fcntl.h>
#include <pthread.h>
#include <stdio.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include <atomic>
#include <fstream>
#include <iostream>
#include <vector>
#define readFileName "/data/test_data.txt"
#define writeFileName "/projects/student/result.txt"
//#define readFileName "./test_data0.txt"
//#define writeFileName "./myresult.txt"
#define N_MAX 4000000
#define E_MAX 2000000
#define buffersize (4096 * 32)
#define bfs_cache_size 40000 * sizeof(bfs3_eum)
#define bfs_cache_sorted_size 40000 * sizeof(bfs3_eum_sorted)
#define THREAD_NUM 4
#define PREVMAX 429496729
#define NEXTMAX 715827882
#include <algorithm>
#include <unordered_map>
using namespace std;
using namespace __gnu_cxx;

typedef struct {
  int nextptr;
  int money_2;
  int id_2;
  int id_1;
  int money_1;
} bfs3_eum;
typedef struct {
  int money_2;
  int id_2;
  int id_1;
  int money_1;
} bfs3_eum_sorted;

typedef struct {
  int threadID;
  int* RadixAddr;
  int* RadixAddrRadix;
  bfs3_eum* bfs_cache;
  bfs3_eum_sorted* bfs_cache_sorted;
  int resultSize;
  int* result3;
  int* result4;
  int* result5;
  int* result6;
  int* result7;
  int result3Pos;
  int result4Pos;
  int result5Pos;
  int result6Pos;
  int result7Pos;
  char* result3Str;
  char* result4Str;
  char* result5Str;
  char* result6Str;
  char* result7Str;
  int* result3Radix;
  int* result4Radix;
  int* result5Radix;
  int* result6Radix;
  int* result7Radix;
} thread_data;
thread_data td[THREAD_NUM];

char buffer[buffersize];
int bufferpos = 0;
int ptr = 0;
int file_size(FILE* fp) {
  if (!fp) return -1;
  fseek(fp, 0L, SEEK_END);
  int size = ftell(fp);
  return size;
}

typedef struct {
  int left;
  int right;
} node_ptr;

typedef struct {
  unsigned int be_id;
  unsigned int to_id;
  unsigned int money;
} bill;
bill bills[E_MAX];
bill tmp_bills[E_MAX];
int bills_ptr = 0;
char* foutbuffer;
int* result3Seq;
int* result4Seq;
int* result5Seq;
int* result6Seq;
int* result7Seq;
int* PrevInts;
int* NextInts;
int PrevIntsPtr = 0;
int NextIntsPtr = 0;
node_ptr* Prev;
int* NextBuildCount;
int PrevIntsRadix[N_MAX] = {0};
int NextIntsRadix[N_MAX] = {0};
unordered_map<unsigned int, int> mymap;
unsigned int sortedNodes[N_MAX];
int sortedNodesPtr = 1;

char* fileBytes;
int fileSize;
void loadFromTxt() {
  char c = -1;
  unsigned int i[3];
  char* ptr = fileBytes;
  char* ed = ptr + fileSize;
  while (ptr < ed) {
    memset(i, 0, 12);
    while (*ptr != ',') {
      i[0] = i[0] * 10 + *ptr++ - '0';
    }
    ptr++;
    while (*ptr != ',') {
      i[1] = i[1] * 10 + *ptr++ - '0';
    }
    ptr++;
    while (*ptr >= '0') {
      i[2] = i[2] * 10 + *ptr++ - '0';
    }
    while (*ptr++ != '\n')
      ;
    bills[bills_ptr].be_id = i[0];
    bills[bills_ptr].to_id = i[1];
    bills[bills_ptr].money = i[2];
    bills_ptr++;
    //		cout << i[0] << "," << i[1] << "," << i[2] << endl;
  }
}
void initDataFromFile(FILE* fp) {
  fileSize = file_size(fp);  //»ñÈ¡ÎÄ¼þ´óÐ¡fileSize
  int fileBytesSize = fileSize / buffersize * buffersize + buffersize +
                      1;  // fileBytesSize×îÐ¡ÖµÎªbufferSize
  cout << readFileName << endl;
  cout << "fileBytesSize: " << fileBytesSize << endl;
  fileBytes = (char*)malloc(fileBytesSize);  //³õÊ¼»¯ÄÚ´æ¿é£¬´óÐ¡ÎªfileBytesSize

  int tmpInt = 0;
  unsigned int i = 0;
  tmpInt = 0;
  fseek(fp, 0, SEEK_SET);
  while (!feof(fp)) {
    i = fread(fileBytes + tmpInt, 1, buffersize, fp);
    tmpInt += i;
  }
  fileBytes[tmpInt] = '\0';
  fileSize = tmpInt;
  loadFromTxt();
}
char* nodeStr[N_MAX] = {0};
char nodeChar[N_MAX][11];

int my_itoa(int num, char* str) {
  char string[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";

  char* ptr = str;
  int i;
  int j;
  int len = 0;
  if (num < 10) {
    *ptr++ = string[num];
    *ptr = '\0';
    return 1;
  }
  while (num) {
    len++;
    *ptr++ = string[num % 10];
    num /= 10;
    if (num < 10) {
      len++;
      *ptr++ = string[num];
      *ptr = '\0';
      break;
    }
  }

  j = ptr - str - 1;

  for (i = 0; i < (ptr - str) / 2; i++) {
    int temp = str[i];
    str[i] = str[j];
    str[j--] = temp;
  }

  return len;
}
void Swap2(int* arr, int low, int high) {
  int temp[2];
  temp[0] = arr[low * 2], temp[1] = arr[low * 2 + 1];
  arr[low * 2] = arr[high * 2], arr[low * 2 + 1] = arr[high * 2 + 1];
  arr[high * 2] = temp[0], arr[high * 2 + 1] = temp[1];
}
int Partition2(int* arr, int low, int high) {
  int base = arr[low * 2];
  while (low < high) {
    while (low < high && arr[high * 2] >= base) {
      high--;
    }
    Swap2(arr, low, high);
    while (low < high && arr[low * 2] <= base) {
      low++;
    }
    Swap2(arr, low, high);
  }
  return low;
}
void QuickSort2(int arr[], int low, int high) {
  if (low < high) {
    int base = Partition2(arr, low, high);
    QuickSort2(arr, low, base - 1);
    QuickSort2(arr, base + 1, high);
  }
}

pthread_barrier_t barrier;
inline bool cmp_by_be_id(const bill& a, const bill& b) {
  if (a.be_id > b.be_id) {
    return false;
  } else if (a.be_id < b.be_id) {
    return true;
  }
  if (a.to_id > b.to_id) {
    return false;
  } else if (a.to_id < b.to_id) {
    return true;
  }
}
inline bool cmp_by_to_id(const bill& a, const bill& b) {
  if (a.to_id >= b.to_id) {
    return false;
  } else {
    return true;
  }
}
void* thread_sort(void* arg) {
  int* task = (int*)arg;
  int beg = task[0];
  int end = task[1];
  sort(bills + beg, bills + end, cmp_by_to_id);
  cout << bills_ptr << " === "
       << " : " << beg << "->" << end << endl;
  pthread_barrier_wait(&barrier);
  pthread_exit(NULL);
}

void meger(int* task_arr) {
  int index[THREAD_NUM];
  memcpy(index, task_arr, THREAD_NUM * sizeof(int));
  bills[bills_ptr].be_id = 1 << 31;
  bills[bills_ptr].to_id = 1 << 31;

  int tmp_to_id = -1;
  for (int i = 0; i < bills_ptr; ++i) {
    int min_index;
    bill min_bill = bills[bills_ptr];
    for (int j = 0; j < THREAD_NUM; ++j) {
      if ((index[j] < task_arr[j + 1]) &&
          cmp_by_to_id(bills[index[j]], min_bill)) {
        min_index = j;
        min_bill = bills[index[j]];
      }
    }
    //        cout << min_bill.be_id << "," << min_bill.to_id << "," <<
    //        min_bill.money << endl;
    if (min_bill.to_id != tmp_to_id) {
      tmp_to_id = min_bill.to_id;
      mymap[tmp_to_id] = sortedNodesPtr;
      sortedNodes[sortedNodesPtr++] = tmp_to_id;
    }
    index[min_index]++;
    tmp_bills[i].be_id = min_bill.be_id;
    tmp_bills[i].to_id = sortedNodesPtr - 1;
    tmp_bills[i].money = min_bill.money;
  }
  //	Prev[tmp_to_id].right = PrevIntsPtr;
}
void* thread_sort2(void* arg) {
  int* task = (int*)arg;
  int beg = task[0];
  int end = task[1];
  sort(tmp_bills + beg, tmp_bills + end, cmp_by_be_id);
  cout << bills_ptr << " === "
       << " : " << beg << "->" << end << endl;
  pthread_barrier_wait(&barrier);
  pthread_exit(NULL);
}
void* build_next_graph(void* arg) {
  int tid = (int)(long long)arg;
  for (int i = 0; i < bills_ptr; i++) {
    bill* tmp_bill = &tmp_bills[i];
    int b0 = tmp_bill->be_id;
    if (b0 % THREAD_NUM != tid || b0 == 0) {
      continue;
    }
    if (NextBuildCount[b0] == 0) {
      char* string = nodeChar[b0];
      my_itoa(sortedNodes[b0], string);
      nodeStr[b0] = string;
      NextBuildCount[b0] = PrevIntsRadix[b0 - 1];
    }
    NextInts[NextBuildCount[b0]++] = tmp_bill->to_id;
    NextInts[NextBuildCount[b0]++] = tmp_bill->money;
  }
  pthread_barrier_wait(&barrier);
  pthread_exit(NULL);
}

void init2() {
  mymap.reserve(bills_ptr);
  pthread_t threads[THREAD_NUM];
  int task_arr[THREAD_NUM + 1] = {0};
  int task_num = bills_ptr / THREAD_NUM;
  for (int i = 0; i < THREAD_NUM; i++) {
    task_arr[i + 1] = task_arr[i] + task_num;
  }
  task_arr[THREAD_NUM] = bills_ptr;
  pthread_barrier_init(&barrier, NULL, THREAD_NUM + 1);
  for (int i = 0; i < THREAD_NUM; i++) {
    int ret = pthread_create(threads + i, NULL, &thread_sort,
                             task_arr + i);  //³É¹¦·µ»Ø0,´íÎó·µ»Ø´íÎó±àºÅ
    if (ret != 0) {
      printf("Create pthread error!\n");
      exit(1);
    }
  }
  pthread_barrier_wait(&barrier);
  meger(task_arr);

  sortedNodesPtr++;
  PrevInts = (int*)malloc(bills_ptr * sizeof(int) * 2);
  NextInts = (int*)malloc(bills_ptr * sizeof(int) * 2);
  result3Seq = (int*)malloc((sortedNodesPtr + 1) * sizeof(int));
  result4Seq = (int*)malloc((sortedNodesPtr + 1) * sizeof(int));
  result5Seq = (int*)malloc((sortedNodesPtr + 1) * sizeof(int));
  result6Seq = (int*)malloc((sortedNodesPtr + 1) * sizeof(int));
  result7Seq = (int*)malloc((sortedNodesPtr + 2) * sizeof(int));
  Prev = (node_ptr*)malloc((sortedNodesPtr + 1) * sizeof(node_ptr));
  NextBuildCount = (int*)malloc((sortedNodesPtr + 1) * sizeof(node_ptr));
  memset(Prev, 0, (sortedNodesPtr + 1) * sizeof(node_ptr));
  memset(NextBuildCount, 0, (sortedNodesPtr + 1) * sizeof(int));
  memset(result3Seq, 0, (sortedNodesPtr + 1) * sizeof(int));
  memset(result4Seq, 0, (sortedNodesPtr + 1) * sizeof(int));
  memset(result5Seq, 0, (sortedNodesPtr + 1) * sizeof(int));
  memset(result6Seq, 0, (sortedNodesPtr + 1) * sizeof(int));
  memset(result7Seq, 0, (sortedNodesPtr + 2) * sizeof(int));

  int tmp_to_id = -1;
  for (int i = 0; i < bills_ptr; i++) {
    bill* tmp_bill = &tmp_bills[i];
    unordered_map<unsigned int, int>::iterator tmp_iter =
        mymap.find(tmp_bill->be_id);
    if (tmp_iter != mymap.end()) {
      tmp_bill->be_id = tmp_iter->second;
      PrevIntsRadix[tmp_bill->be_id] += 2;
      if (tmp_to_id != tmp_bill->to_id) {
        if (tmp_to_id != -1) {
          Prev[tmp_to_id].right = PrevIntsPtr;
        }
        tmp_to_id = tmp_bill->to_id;
        Prev[tmp_to_id].left = PrevIntsPtr;
      }
      PrevInts[PrevIntsPtr++] = tmp_bill->be_id;
      PrevInts[PrevIntsPtr++] = tmp_bill->money;
    } else {
      tmp_bill->be_id = 0;
    }
  }
  if (tmp_to_id != -1) {
    NextIntsRadix[tmp_to_id] = PrevIntsPtr;
    Prev[tmp_to_id].right = PrevIntsPtr;
  }

  for (int i = 0; i < sortedNodesPtr - 1; i++) {
    PrevIntsRadix[i + 1] += PrevIntsRadix[i];
    if (Prev[i + 1].right == 0) {
      NextIntsRadix[i + 1] = NextIntsRadix[i];
      continue;
    }
    NextIntsRadix[i] = Prev[i + 1].left;
    NextIntsRadix[i + 1] = Prev[i + 1].right;
  }

  pthread_barrier_init(&barrier, NULL, THREAD_NUM + 1);
  for (int i = 0; i < THREAD_NUM; i++) {
    int ret = pthread_create(threads + i, NULL, &build_next_graph,
                             (void*)(i));  //³É¹¦·µ»Ø0,´íÎó·µ»Ø´íÎó±àºÅ
    if (ret != 0) {
      printf("Create pthread error!\n");
      exit(1);
    }
  }
  pthread_barrier_wait(&barrier);
}
inline bool cmp(int X, int Y) { return X <= 5ul * Y && Y <= 3ul * X; }
bool cmp_bfs3(const bfs3_eum_sorted& a, const bfs3_eum_sorted& b) {
  if (a.id_2 > b.id_2) {
    return false;
  } else if (a.id_2 < b.id_2) {
    return true;
  }
  if (a.id_1 > b.id_1) {
    return false;
  } else if (a.id_1 < b.id_1) {
    return true;
  }
}
inline int sort_bfs3(bfs3_eum* bfs_cache, bfs3_eum_sorted* bfs_cache_sorted,
                     int& bfs_cache_ptr, int ptr) {
  int tmp_ptr = bfs_cache_ptr++;
  bfs_cache_sorted[tmp_ptr].id_2 = 0;
  do {
    bfs_cache_sorted[bfs_cache_ptr].money_2 = bfs_cache[ptr].money_2;
    bfs_cache_sorted[bfs_cache_ptr].id_2 = bfs_cache[ptr].id_2;
    bfs_cache_sorted[bfs_cache_ptr].id_1 = bfs_cache[ptr].id_1;
    bfs_cache_sorted[bfs_cache_ptr].money_1 = bfs_cache[ptr].money_1;
    bfs_cache_ptr++;
    ptr = bfs_cache[ptr].nextptr;
  } while (ptr != 0);
  sort(&bfs_cache_sorted[tmp_ptr], &bfs_cache_sorted[bfs_cache_ptr], cmp_bfs3);
  bfs_cache_sorted[tmp_ptr].money_2 = bfs_cache_ptr;
  return -tmp_ptr;
}
inline char* sync_char(char* tmpChar, int point) {
  char* strSrc = nodeStr[point];
  while ((*tmpChar++ = *strSrc++) != '\0')
    ;
  *(tmpChar - 1) = ',';
  return tmpChar;
}
int getResult(int htoT, thread_data* td, int resultStr[]) {
  if (NextBuildCount[htoT] == 0) {
    return 0;
  }
  int resultSize = 0;
  int* tmp;

  bfs3_eum* bfs_cache = td->bfs_cache;
  int bfs_cache_ptr = 1;
  int* RadixAddr = td->RadixAddr;
  int* RadixAddrRadix = td->RadixAddrRadix;
  int RadixAddrRadixPos = 0;

  int pass = 1;
  register int i0, i1, i2, i3, i4, i5, i6, i7;

  int point[7];
  point[0] = htoT;

  int tmpIntMoney[7];

  for (i0 = NextIntsRadix[point[0] - 1]; i0 < NextIntsRadix[point[0]];
       i0 += 2) {
    point[1] = PrevInts[i0];
    if (point[1] <= point[0]) {
      continue;
    }
    tmpIntMoney[0] = PrevInts[i0 + 1];
    for (i1 = NextIntsRadix[point[1] - 1]; i1 < NextIntsRadix[point[1]];
         i1 += 2) {
      point[2] = PrevInts[i1];
      if (point[2] <= point[0]) {
        continue;
      }
      tmpIntMoney[1] = PrevInts[i1 + 1];

      if (cmp(tmpIntMoney[1], tmpIntMoney[0]) == false) {
        continue;
      }

      for (i2 = NextIntsRadix[point[2] - 1]; i2 < NextIntsRadix[point[2]];
           i2 += 2) {
        point[3] = PrevInts[i2];
        if (point[3] < point[0] || point[3] == point[1]) {
          continue;
        }
        pass = 0;
        if (point[3] == point[0]) {
          continue;
        }
        tmpIntMoney[2] = PrevInts[i2 + 1];
        if (cmp(tmpIntMoney[2], tmpIntMoney[1]) == false) {
          continue;
        }
        if (RadixAddr[point[3]] == 0) {
          RadixAddr[point[3]] = bfs_cache_ptr;
          bfs_cache[bfs_cache_ptr].nextptr = 0;
          bfs_cache[bfs_cache_ptr].money_2 = tmpIntMoney[2];
          bfs_cache[bfs_cache_ptr].id_2 = point[2];
          bfs_cache[bfs_cache_ptr].id_1 = point[1];
          bfs_cache[bfs_cache_ptr].money_1 = tmpIntMoney[0];
          bfs_cache_ptr++;
          RadixAddrRadix[RadixAddrRadixPos++] = point[3];
        } else {
          int tmp_RadixAddr = RadixAddr[point[3]];
          RadixAddr[point[3]] = bfs_cache_ptr;
          bfs_cache[bfs_cache_ptr].nextptr = tmp_RadixAddr;
          bfs_cache[bfs_cache_ptr].money_2 = tmpIntMoney[2];
          bfs_cache[bfs_cache_ptr].id_2 = point[2];
          bfs_cache[bfs_cache_ptr].id_1 = point[1];
          bfs_cache[bfs_cache_ptr].money_1 = tmpIntMoney[0];
          bfs_cache_ptr++;
        }
      }
    }
  }
  if (pass) {
    int i = 0;
    while (i < RadixAddrRadixPos) {
      RadixAddr[RadixAddrRadix[i]] = 0;
      i++;
    }
    return 0;
  }

  int* result3 = td->result3;
  int* result4 = td->result4;
  int* result5 = td->result5;
  int* result6 = td->result6;
  int* result7 = td->result7;
  char* result3Str = td->result3Str;
  char* result4Str = td->result4Str;
  char* result5Str = td->result5Str;
  char* result6Str = td->result6Str;
  char* result7Str = td->result7Str;
  int result3Pos = td->result3Pos;
  int result4Pos = td->result4Pos;
  int result5Pos = td->result5Pos;
  int result6Pos = td->result6Pos;
  int result7Pos = td->result7Pos;
  int result3StrPos = resultStr[0];
  int result4StrPos = resultStr[1];
  int result5StrPos = resultStr[2];
  int result6StrPos = resultStr[3];
  int result7StrPos = resultStr[4];
  int* result3Radix = td->result3Radix;
  int* result4Radix = td->result4Radix;
  int* result5Radix = td->result5Radix;
  int* result6Radix = td->result6Radix;
  int* result7Radix = td->result7Radix;

  bfs_cache_ptr = 1;
  bfs3_eum_sorted* bfs_cache_sorted = td->bfs_cache_sorted;

  int int_to_char_result[7] = {-1};
  char tmpCharResult[7 * 11];
  char* tmpChar[8];
  tmpChar[0] = tmpCharResult;

  for (i0 = PrevIntsRadix[point[0] - 1]; i0 < PrevIntsRadix[point[0]];
       i0 += 2) {
    point[1] = NextInts[i0];
    if (point[1] < point[0]) {
      continue;
    }

    tmpIntMoney[0] = NextInts[i0 + 1];

    for (i1 = PrevIntsRadix[point[1] - 1]; i1 < PrevIntsRadix[point[1]];
         i1 += 2) {
      point[2] = NextInts[i1];
      if (point[2] <= point[0]) {
        continue;
      }
      tmpIntMoney[1] = NextInts[i1 + 1];
      if (cmp(tmpIntMoney[0], tmpIntMoney[1]) == false) {
        continue;
      }
      //²éÎå½Úµã»·
      if (RadixAddr[point[2]] != 0) {
        int ptr = RadixAddr[point[2]], ptr2;
        if (ptr > 0) {
          ptr = sort_bfs3(bfs_cache, bfs_cache_sorted, bfs_cache_ptr, ptr);
          RadixAddr[point[2]] = ptr;
        }
        ptr = -ptr;
        ptr2 = bfs_cache_sorted[ptr++].money_2;
        for (; ptr < ptr2; ptr++) {
          if (cmp(tmpIntMoney[1], bfs_cache_sorted[ptr].money_2) == false ||
              cmp(bfs_cache_sorted[ptr].money_1, tmpIntMoney[0]) == false) {
            continue;
          }
          point[3] = bfs_cache_sorted[ptr].id_2;
          point[4] = bfs_cache_sorted[ptr].id_1;

          int p = 1, q = 0, pq = 0;
          for (; p < 3; p++) {
            for (q = 3; q < 5; q++) {
              if (point[p] == point[q]) {
                pq++;
              }
            }
          }
          if (pq) {
            continue;
          }
          for (int tmpi = 0; tmpi < 5; tmpi++) {
            if (int_to_char_result[tmpi] != point[tmpi]) {
              for (int i = tmpi; i < 5; i++) {
                int_to_char_result[i] = point[i];
                tmpChar[i + 1] = tmpChar[i];
                tmpChar[i + 1] = sync_char(tmpChar[i + 1], point[i]);
              }
              int_to_char_result[5] = -1;
              break;
            }
          }
          int strlen = tmpChar[5] - tmpCharResult;
          memcpy(result5Str + result5StrPos, tmpCharResult, strlen);
          result5StrPos += strlen;
          *(result5Str + result5StrPos - 1) = '\n';
          resultSize++;
        }
      }
      for (i2 = PrevIntsRadix[point[2] - 1]; i2 < PrevIntsRadix[point[2]];
           i2 += 2) {
        point[3] = NextInts[i2];
        if (point[3] < point[0] || point[3] == point[1]) {
          continue;
        }
        tmpIntMoney[2] = NextInts[i2 + 1];
        if (cmp(tmpIntMoney[1], tmpIntMoney[2]) == false) {
          continue;
        }

        //²éÈý½Úµã»·
        if (point[3] == point[0]) {
          if (cmp(tmpIntMoney[2], tmpIntMoney[0])) {
            for (int tmpi = 0; tmpi < 3; tmpi++) {
              if (int_to_char_result[tmpi] != point[tmpi]) {
                for (int i = tmpi; i < 3; i++) {
                  int_to_char_result[i] = point[i];
                  tmpChar[i + 1] = tmpChar[i];
                  tmpChar[i + 1] = sync_char(tmpChar[i + 1], point[i]);
                }
                int_to_char_result[3] = -1;
                break;
              }
            }
            int strlen = tmpChar[3] - tmpCharResult;
            memcpy(result3Str + result3StrPos, tmpCharResult, strlen);
            result3StrPos += strlen;
            *(result3Str + result3StrPos - 1) = '\n';
            resultSize++;
          }
          continue;
        }
        //²éÁù½Úµã»·
        if (RadixAddr[point[3]] != 0) {
          int ptr = RadixAddr[point[3]], ptr2;
          if (ptr > 0) {
            ptr = sort_bfs3(bfs_cache, bfs_cache_sorted, bfs_cache_ptr, ptr);
            RadixAddr[point[3]] = ptr;
          }
          ptr = -ptr;
          ptr2 = bfs_cache_sorted[ptr++].money_2;
          for (; ptr < ptr2; ptr++) {
            if (cmp(tmpIntMoney[2], bfs_cache_sorted[ptr].money_2) == false ||
                cmp(bfs_cache_sorted[ptr].money_1, tmpIntMoney[0]) == false) {
              continue;
            }
            point[4] = bfs_cache_sorted[ptr].id_2;
            point[5] = bfs_cache_sorted[ptr].id_1;

            int p = 1, q = 0, pq = 0;
            for (; p < 4; p++) {
              for (q = 4; q < 6; q++) {
                if (point[p] == point[q]) {
                  pq++;
                }
              }
            }
            if (pq) {
              continue;
            }
            for (int tmpi = 0; tmpi < 6; tmpi++) {
              if (int_to_char_result[tmpi] != point[tmpi]) {
                for (int i = tmpi; i < 6; i++) {
                  int_to_char_result[i] = point[i];
                  tmpChar[i + 1] = tmpChar[i];
                  tmpChar[i + 1] = sync_char(tmpChar[i + 1], point[i]);
                }
                int_to_char_result[6] = -1;
                break;
              }
            }
            int strlen = tmpChar[6] - tmpCharResult;
            memcpy(result6Str + result6StrPos, tmpCharResult, strlen);
            result6StrPos += strlen;
            *(result6Str + result6StrPos - 1) = '\n';
            resultSize++;
          }
        }
        //				exit(1);
        for (i3 = PrevIntsRadix[point[3] - 1]; i3 < PrevIntsRadix[point[3]];
             i3 += 2) {
          point[4] = NextInts[i3];
          tmpIntMoney[3] = NextInts[i3 + 1];
          //²éËÄ½Úµã»·
          if (point[4] == point[0]) {
            if (point[4] == point[1] || point[4] == point[2]) {
              continue;
            }
            if (cmp(tmpIntMoney[2], tmpIntMoney[3]) == false) {
              continue;
            }
            if (cmp(tmpIntMoney[3], tmpIntMoney[0])) {
              for (int tmpi = 0; tmpi < 4; tmpi++) {
                if (int_to_char_result[tmpi] != point[tmpi]) {
                  for (int i = tmpi; i < 4; i++) {
                    int_to_char_result[i] = point[i];
                    tmpChar[i + 1] = tmpChar[i];
                    tmpChar[i + 1] = sync_char(tmpChar[i + 1], point[i]);
                  }
                  int_to_char_result[4] = -1;
                  break;
                }
              }
              int strlen = tmpChar[4] - tmpCharResult;
              memcpy(result4Str + result4StrPos, tmpCharResult, strlen);
              result4StrPos += strlen;
              *(result4Str + result4StrPos - 1) = '\n';
              resultSize++;
            }
            continue;
          } else if (point[4] < point[0]) {
            continue;
          }

          //²éÆß½Úµã»·
          if (RadixAddr[point[4]] != 0) {
            if (point[4] == point[1] || point[4] == point[2]) {
              continue;
            }
            if (cmp(tmpIntMoney[2], tmpIntMoney[3]) == false) {
              continue;
            }
            int ptr = RadixAddr[point[4]], ptr2;
            if (ptr > 0) {
              ptr = sort_bfs3(bfs_cache, bfs_cache_sorted, bfs_cache_ptr, ptr);
              RadixAddr[point[4]] = ptr;
            }
            ptr = -ptr;
            ptr2 = bfs_cache_sorted[ptr++].money_2;
            for (; ptr < ptr2; ptr++) {
              if (cmp(tmpIntMoney[3], bfs_cache_sorted[ptr].money_2) == false ||
                  cmp(bfs_cache_sorted[ptr].money_1, tmpIntMoney[0]) == false) {
                continue;
              }
              point[5] = bfs_cache_sorted[ptr].id_2;
              point[6] = bfs_cache_sorted[ptr].id_1;

              int p = 1, q = 0, pq = 0;
              for (; p < 5; p++) {
                for (q = 5; q < 7; q++) {
                  if (point[p] == point[q]) {
                    pq++;
                  }
                }
              }
              if (pq) {
                continue;
              }

              for (int tmpi = 0; tmpi < 7; tmpi++) {
                if (int_to_char_result[tmpi] != point[tmpi]) {
                  for (int i = tmpi; i < 7; i++) {
                    int_to_char_result[i] = point[i];
                    tmpChar[i + 1] = tmpChar[i];
                    tmpChar[i + 1] = sync_char(tmpChar[i + 1], point[i]);
                  }
                  break;
                }
              }
              int strlen = tmpChar[7] - tmpCharResult;
              memcpy(result7Str + result7StrPos, tmpCharResult, strlen);
              result7StrPos += strlen;
              *(result7Str + result7StrPos - 1) = '\n';
              resultSize++;
            }
          }
        }
      }
    }
  }
  //	cout << "==" << endl;
  if (resultStr[0] != result3StrPos) {
    int Pos = result3Pos * 3;
    int len = result3StrPos - resultStr[0];
    result3Radix[result3Pos] = Pos;
    result3Seq[point[0] + 1] = len;
    result3[Pos] = point[0];
    result3[Pos + 1] = resultStr[0];
    result3[Pos + 2] = len;
    result3Pos++;
    td->result3Pos = result3Pos;
    resultStr[0] = result3StrPos;
  }
  if (resultStr[1] != result4StrPos) {
    int Pos = result4Pos * 3;
    int len = result4StrPos - resultStr[1];
    result4Radix[result4Pos] = Pos;
    result4Seq[point[0] + 1] = len;
    result4[Pos] = point[0];
    result4[Pos + 1] = resultStr[1];
    result4[Pos + 2] = len;
    result4Pos++;
    td->result4Pos = result4Pos;
    resultStr[1] = result4StrPos;
  }
  if (resultStr[2] != result5StrPos) {
    int Pos = result5Pos * 3;
    int len = result5StrPos - resultStr[2];
    result5Radix[result5Pos] = Pos;
    result5Seq[point[0] + 1] = len;
    result5[Pos] = point[0];
    result5[Pos + 1] = resultStr[2];
    result5[Pos + 2] = len;
    result5Pos++;
    td->result5Pos = result5Pos;
    resultStr[2] = result5StrPos;
  }
  if (resultStr[3] != result6StrPos) {
    int Pos = result6Pos * 3;
    int len = result6StrPos - resultStr[3];
    result6Radix[result6Pos] = Pos;
    result6Seq[point[0] + 1] = len;
    result6[Pos] = point[0];
    result6[Pos + 1] = resultStr[3];
    result6[Pos + 2] = len;
    result6Pos++;
    td->result6Pos = result6Pos;
    resultStr[3] = result6StrPos;
  }
  if (resultStr[4] != result7StrPos) {
    int Pos = result7Pos * 3;
    int len = result7StrPos - resultStr[4];
    result7Radix[result7Pos] = Pos;
    result7Seq[point[0] + 1] = len;
    result7[Pos] = point[0];
    result7[Pos + 1] = resultStr[4];
    result7[Pos + 2] = len;
    result7Pos++;
    td->result7Pos = result7Pos;
    resultStr[4] = result7StrPos;
  }
  int i = 0;
  while (i < RadixAddrRadixPos) {
    RadixAddr[RadixAddrRadix[i]] = 0;
    i++;
  }
  return resultSize;
}

atomic<int> atomic_i(1);
void* handleHtoT(void* threadarg) {
  int threadID = (int)(long long)threadarg;
  int resultSize = 0;
  int resultStrPos[5] = {0, 0, 0, 0, 0};
  const int fp = 100;
  printf("threadID: %d\n", threadID);
  int i = atomic_i.fetch_add(fp, std::memory_order_relaxed);
  for (; i < sortedNodesPtr;) {
    for (int o = 0; o < fp; ++o) {
      resultSize += getResult(i + o, &td[threadID - 1], resultStrPos);
    }
    i = atomic_i.fetch_add(fp, std::memory_order_relaxed);
  }
  printf("Finish!\n");
  td[threadID - 1].resultSize = resultSize;
  pthread_exit(NULL);
}
void* writeFileBuffer(void* threadarg) {
  thread_data* tmp = (thread_data*)threadarg;
  int o;
  o = 0;
  while (o < tmp->result3Pos) {
    int radix, len, id;
    unsigned startAddr;
    radix = tmp->result3Radix[o];
    if (radix != -1) {
      id = tmp->result3[radix];
      startAddr = tmp->result3[radix + 1];
      len = tmp->result3[radix + 2];
      memcpy(foutbuffer + result3Seq[id], tmp->result3Str + startAddr, len);
    }
    o++;
  }

  o = 0;
  while (o < tmp->result4Pos) {
    int radix, len, id;
    unsigned startAddr;
    radix = tmp->result4Radix[o];
    if (radix != -1) {
      id = tmp->result4[radix];
      startAddr = tmp->result4[radix + 1];
      len = tmp->result4[radix + 2];
      memcpy(foutbuffer + result4Seq[id], tmp->result4Str + startAddr, len);
    }
    o++;
  }
  o = 0;
  while (o < tmp->result5Pos) {
    int radix, len, id;
    unsigned startAddr;
    radix = tmp->result5Radix[o];
    if (radix != -1) {
      id = tmp->result5[radix];
      startAddr = tmp->result5[radix + 1];
      len = tmp->result5[radix + 2];
      memcpy(foutbuffer + result5Seq[id], tmp->result5Str + startAddr, len);
    }
    o++;
  }
  o = 0;
  while (o < tmp->result6Pos) {
    int radix, len, id;
    unsigned startAddr;
    radix = tmp->result6Radix[o];
    if (radix != -1) {
      id = tmp->result6[radix];
      startAddr = tmp->result6[radix + 1];
      len = tmp->result6[radix + 2];
      memcpy(foutbuffer + result6Seq[id], tmp->result6Str + startAddr, len);
    }
    o++;
  }
  o = 0;
  while (o < tmp->result7Pos) {
    int radix, len, id;
    unsigned startAddr;
    radix = tmp->result7Radix[o];
    if (radix != -1) {
      id = tmp->result7[radix];
      startAddr = tmp->result7[radix + 1];
      len = tmp->result7[radix + 2];
      memcpy(foutbuffer + result7Seq[id], tmp->result7Str + startAddr, len);
    }
    o++;
  }
  pthread_exit(NULL);
}
int fout;
int main() {
  FILE* fp = NULL;
  fp = fopen(readFileName, "r");
  if (fp == NULL) {
    printf("file open failed");
    return 0;
  }
  clock_t start = clock();
  initDataFromFile(fp);
  init2();
  cout << clock() - start << endl;
  //	exit(1);
  fclose(fp);
  cout << "hello world!" << endl;
  cout << mymap.size() << endl;

  for (int i = 0; i < THREAD_NUM; i++) {
    td[i].threadID = i + 1;
    int len1 = sortedNodesPtr * sizeof(int);
    int len2 = sortedNodesPtr * sizeof(int);

    int len4 = (sortedNodesPtr + 1) * 3 * sizeof(int);
    int len5 = len4 * 5;
    void* threadBytes = (void*)malloc(len1 + len2 + bfs_cache_size +
                                      bfs_cache_sorted_size + len5);
    td[i].resultSize = 0;
    td[i].RadixAddr = (int*)(threadBytes);
    td[i].RadixAddrRadix = (int*)(threadBytes + len1);
    td[i].bfs_cache = (bfs3_eum*)(threadBytes + len1 + len2);
    td[i].bfs_cache_sorted =
        (bfs3_eum_sorted*)(threadBytes + len1 + len2 + bfs_cache_size);
    td[i].result3 = (int*)(threadBytes + len1 + len2 + bfs_cache_size +
                           bfs_cache_sorted_size);
    td[i].result4 = (int*)(threadBytes + len1 + len2 + bfs_cache_size +
                           bfs_cache_sorted_size + len4);
    td[i].result5 = (int*)(threadBytes + len1 + len2 + bfs_cache_size +
                           bfs_cache_sorted_size + len4 * 2);
    td[i].result6 = (int*)(threadBytes + len1 + len2 + bfs_cache_size +
                           bfs_cache_sorted_size + len4 * 3);
    td[i].result7 = (int*)(threadBytes + len1 + len2 + bfs_cache_size +
                           bfs_cache_sorted_size + len4 * 4);

    void* threadBytes2 = (void*)malloc(len5);
    td[i].result3Radix = (int*)(threadBytes2);
    td[i].result4Radix = (int*)(threadBytes2 + len4);
    td[i].result5Radix = (int*)(threadBytes2 + len4 * 2);
    td[i].result6Radix = (int*)(threadBytes2 + len4 * 3);
    td[i].result7Radix = (int*)(threadBytes2 + len4 * 4);

    memset(td[i].RadixAddr, 0, len1);

    td[i].result3Str = (char*)malloc(1 << 28);
    td[i].result4Str = (char*)malloc(1 << 28);
    td[i].result5Str = (char*)malloc(1 << 28);
    td[i].result6Str = (char*)malloc(1 << 28);
    td[i].result7Str = (char*)malloc(1 << 29);

    td[i].result3Pos = 0;
    td[i].result4Pos = 0;
    td[i].result5Pos = 0;
    td[i].result6Pos = 0;
    td[i].result7Pos = 0;
  }
  pthread_t threads[THREAD_NUM];
  for (int i = 0; i < THREAD_NUM; i++) {
    td[i].threadID = i + 1;
    int ret = pthread_create(&threads[i], NULL, &handleHtoT,
                             (void*)(i + 1));  //³É¹¦·µ»Ø0,´íÎó·µ»Ø´íÎó±àºÅ
    if (ret != 0) {
      printf("Create pthread error!\n");
      exit(1);
    }
  }
  int allresultSize = 0;
  for (int i = 0; i < THREAD_NUM; i++) {
    pthread_join(threads[i], NULL);
    allresultSize += td[i].resultSize;
  }
  cout << "allresultSize: " << allresultSize << endl;

  printf("Start Write!\n");

  //	-------------------Ð´ÎÄ¼þ
  // exit(1);

  fout = open(writeFileName, O_RDWR | O_CREAT, 0666);
  if (fout == NULL) {
    printf("file open failed");
    exit(1);
  }

  char aaa[12];
  int o = 1;
  result3Seq[0] = my_itoa(allresultSize, aaa) + 1;
  while (o <= sortedNodesPtr) {
    result3Seq[o] += result3Seq[o - 1];
    o++;
  }
  result4Seq[0] += result3Seq[o - 1];
  o = 1;
  while (o <= sortedNodesPtr) {
    result4Seq[o] += result4Seq[o - 1];
    o++;
  }
  result5Seq[0] += result4Seq[o - 1];
  o = 1;
  while (o <= sortedNodesPtr) {
    result5Seq[o] += result5Seq[o - 1];
    o++;
  }
  result6Seq[0] += result5Seq[o - 1];
  o = 1;
  while (o <= sortedNodesPtr) {
    result6Seq[o] += result6Seq[o - 1];
    o++;
  }
  result7Seq[0] += result6Seq[o - 1];
  o = 1;
  while (o <= sortedNodesPtr + 1) {
    result7Seq[o] += result7Seq[o - 1];
    o++;
  }
  int allfilesize = result7Seq[sortedNodesPtr + 1];
  cout << "allfilesize: " << allfilesize << endl;

  truncate(writeFileName, allfilesize);
  foutbuffer = (char*)mmap(NULL, allfilesize, PROT_READ | PROT_WRITE,
                           MAP_SHARED, fout, 0);
  allresultSize = my_itoa(allresultSize, foutbuffer);
  *(foutbuffer + allresultSize) = '\n';
  //	---------------------------Ë¢ÐÂÎÄ¼þ
  // exit(1);

  for (int i = 0; i < THREAD_NUM; i++) {
    int ret = pthread_create(threads + i, NULL, &writeFileBuffer,
                             td + i);  //³É¹¦·µ»Ø0,´íÎó·µ»Ø´íÎó±àºÅ
    if (ret != 0) {
      printf("Create pthread error!\n");
      while (1) {
      }
    }
  }
  for (int i = 0; i < THREAD_NUM; i++) {
    pthread_join(threads[i], NULL);
  }
  munmap(foutbuffer, allfilesize);
}
