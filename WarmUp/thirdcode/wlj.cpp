#include <fcntl.h>     /* for open */
#include <sys/mman.h>  /* for mmap and munmap */
#include <sys/stat.h>  /* for open */
#include <sys/types.h> /* for open */
#include <unistd.h>    /* for lseek and write */
#include <cmath>
#include <condition_variable>
#include <ctime>
#include <fstream>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

using namespace std;
mutex g_mutex_weight;
mutex g_mutex_predict;
condition_variable cv_weight;
bool isWeight = false;
condition_variable cv_predict;
const int predictNum = 3;
vector<int> predictPos(predictNum + 1, -1);
vector<mutex> g_mutex_predictPos(predictNum + 1);
vector<condition_variable> cv_predictPos(predictNum + 1);

typedef float dataType;

inline void sigmoid(dataType &value) { value = 1 / (1 + exp(-value)); }

void nomalization(vector<vector<dataType>> &data_train,
                  vector<vector<dataType>>::iterator weight,
                  int arrayNum)  //标准化处理
{
  for (int i = 0; i < arrayNum; i++) {
    for (int j = 0; j < 1000; j++)

      data_train[i][j] = (data_train[i][j] - weight[1][j]) / sqrt(weight[2][j]);
  }
}

inline void opV(vector<dataType> &vector_v, const vector<dataType> &vector_c,
                const dataType &coe_v = 1, const dataType &coe_c = 1,
                int index = 1000)  // train函数的子函数
{
  while (index-- > 0) vector_v[index] += vector_c[index] * coe_c;
}

void mul_row_col(const vector<dataType> &row, const vector<dataType> &col,
                 vector<dataType> &grad_weight, int batch)  // train函数的子函数
{
  int index = 1000;
  dataType result = 0;
  while (index-- > 0) {
    result += row[index] * col[index];
    /* code */
  }
  index = 1002;
  result += col[index - 1];
  sigmoid(result);

  opV(grad_weight, row, 1, (result - row[index - 2]) / batch, index);
}

void train(const vector<vector<dataType>> &data_train, vector<dataType> &weight,
           int data_num)  // lr训练数据
{
  int batchsize = 32;
  int batchEpoch = data_num / batchsize;
  int maxEpoch = 2;
  int epoch = 0;
  vector<dataType> grad_weight(1024);
  double learningRate = -0.4;

  while (epoch++ < maxEpoch) {
    int index = -1;
    while (++index < batchEpoch) {
      int batch = -1;
      int dim = 1000;
      while (++batch < batchsize) {
        int inputIndex = index * batchsize + batch;
        mul_row_col(data_train[inputIndex], weight, grad_weight, batchsize);
      }
      // opV(grad_weight, weight, 1, lamdb / batchsize);
      opV(weight, grad_weight, 0, learningRate, dim + 2);
    }
    if (epoch == 0) learningRate = -0.06;
  }
}

void readFile(const char *filePath,
              vector<vector<dataType>>::iterator
                  weight)  //读训练数据，读完后调用train函数训练
{
  // clock_t start_time = clock();
  const int arrayNumS = 1900;
  const int fileSize = arrayNumS * 1000 * 6.5;
  int arrayNum = arrayNumS;
  vector<vector<dataType>> x_input(arrayNum + 1, vector<dataType>(1024));

  // C++标准库读文件
  // char *const buf = new char[fileSize];
  // ifstream fin(filePath, ios::binary);
  // fin.read(&buf[0], static_cast<streamsize>(fileSize));
  // linux读文件
  int fd = open(filePath, O_RDONLY);
  char *buf = (char *)mmap(NULL, fileSize, PROT_READ, MAP_PRIVATE, fd, 0);

  char *iter_fileStart = buf;
  auto iter_row = x_input.begin();
  vector<dataType>::iterator iter_col = x_input[0].begin();
  const vector<dataType>::iterator iter_mean_start = weight[1].begin();
  vector<dataType>::iterator iter_mean = iter_mean_start;
  const vector<dataType>::iterator iter_var_start = weight[2].begin();
  vector<dataType>::iterator iter_var = iter_var_start;
  while (arrayNum > 0) {
    bool neg = false;
    if (*iter_fileStart == '-') {
      neg = true;
      ++iter_fileStart;
    }
    double f = *iter_fileStart++ - '0';

    if (*iter_fileStart == '.') {
      ++iter_fileStart;
      while (*iter_fileStart >= '0' && *iter_fileStart <= '9') {
        f = (f * 10.0) + (*iter_fileStart - '0');
        ++iter_fileStart;
        // ++n;
      }
      f = f * 0.001;
    }

    if (neg) {
      f = -f;
    }
    *iter_col++ = f;
    *iter_mean++ += f / arrayNumS;
    if (*iter_fileStart == ',') {
      /* code */
      iter_fileStart++;
      // continue;
    } else if (*iter_fileStart == '\n') {
      iter_mean = iter_mean_start;
      // arrayNum++;
      *(iter_col) = 1;
      iter_fileStart++;
      iter_row++;
      iter_col = (*iter_row).begin();
      arrayNum--;
      // if (arrayNum ==0)
      // break;

      // continue;
    }
  }

  // delete buf;
  close(fd);
  munmap(buf, fileSize);
  // clock_t end_time = clock();

  iter_row = x_input.begin();
  iter_col = (*iter_row).begin();
  arrayNum = arrayNumS;
  int dim = 1000;
  double temp = 0;
  while (arrayNum > 0) {
    temp = *iter_mean - *iter_col;
    iter_col++;
    iter_mean++;
    *iter_var++ += temp * temp / (arrayNumS - 1);
    dim--;
    if (!dim) {
      dim = 1000;
      iter_mean = iter_mean_start;
      iter_row++;
      iter_col = (*iter_row).begin();
      iter_var = iter_var_start;
      arrayNum--;
    }
  }
  nomalization(x_input, weight, arrayNumS);

  // cout << "\nloadData time is: " << (double)(end_time - start_time) /
  // CLOCKS_PER_SEC << "s" << endl;

  // start_time = clock();
  train(x_input, weight[0], arrayNumS);
  // end_time = clock();
  // cout << "\nThe train time is: " << (double)(end_time - start_time) /
  // CLOCKS_PER_SEC << "s" << endl;
}

inline char mulV(const vector<dataType> &row,
                 vector<vector<dataType>>::iterator weight)  // precdict的子函数
{
  int index = 1000;
  dataType result = 0;
  while (index-- > 0)
    result += (row[index] - weight[1][index]) * weight[0][index] /
              sqrt(weight[2][index]);

  return (result + weight[0][1001]) < 0 ? '0' : '1';
}

void predict(vector<char>::iterator predictData,
             const vector<vector<dataType>> &data,
             vector<vector<dataType>>::iterator weight, int data_num,
             int threadId)  //预测测试数据
{
  // cout << predictPos[1] << endl;
  // const int bufferSize = 2 * data_num;
  int start = predictPos[threadId];
  int index = -1;
  while (++index < data_num)
    predictData[2 * (start + index)] = mulV(data[index], weight);
  // cout << endl;
  // cout<<"write predict
  // testData:"<<(clock()-start_time)/(doubl=e)CLOCKS_PER_SEC;
}

void analData_thread(vector<char>::iterator predictData, char *iter_fileStart,
                     char *iter_fileEnd,
                     vector<vector<dataType>>::iterator weight,
                     int threadId)  //解析测试数据
{
  if (threadId) {
    if (*(iter_fileStart - 1) != '\n')
      while (*iter_fileStart++ != '\n')
        ;
  }
  vector<vector<dataType>> x_input((iter_fileEnd - iter_fileStart) / 6000 + 1,
                                   vector<dataType>(1024));
  auto iter_row = x_input.begin();
  auto iter_col = (*iter_row).begin();
  int arrayNum = 0;
  while (true) {
    bool neg = false;
    if (*iter_fileStart == '-') {
      neg = true;
      ++iter_fileStart;
    }
    double f = *iter_fileStart++ - '0';
    // if (*iter_fileStart == '.')
    // {
    /* code */
    // double f = r;
    // int n = 0;
    ++iter_fileStart;
    while (*iter_fileStart >= '0' && *iter_fileStart <= '9') {
      f = (f * 10.0) + (*iter_fileStart - '0');
      ++iter_fileStart;
      // ++n;
    }
    f = f * 0.001;  // pow(0.1, 3);
    // }

    if (neg) {
      f = -f;
    }
    *iter_col = f;
    if (*iter_fileStart == ',') {
      /* code */
      iter_fileStart++;
      iter_col++;
      // continue;
    } else if (*iter_fileStart == '\n') {
      arrayNum++;
      iter_fileStart++;
      if (iter_fileStart >= iter_fileEnd) break;
      iter_row++;
      iter_col = (*iter_row).begin();
    }
  }
  // cout << "point1" << endl;
  unique_lock<mutex> locker_weight(g_mutex_weight);
  if (!isWeight) {
    cv_weight.wait(locker_weight);
    locker_weight.unlock();
    cv_weight.notify_all();
  } else
    locker_weight.unlock();
  unique_lock<mutex> locker_pos_self(g_mutex_predictPos[threadId]);
  if (predictPos[threadId] < 0) {
    cv_predictPos[threadId].wait(locker_pos_self);
  }
  {
    unique_lock<mutex> locker_pos_next(g_mutex_predictPos[threadId + 1]);
    predictPos[threadId + 1] = predictPos[threadId] + arrayNum;
  }
  cv_predictPos[threadId + 1].notify_one();

  predict(predictData, x_input, weight, arrayNum, threadId);
}

void readtestFile(
    const char *filePath, const char *filePath_save,
    vector<vector<dataType>>::iterator weight)  //读测试数据，读完后预测
{
  // unique_lock<mutex> locker1(g_mutex_predict);
  // const char* filePath = "data/test_data_m.txt";
  // clock_t start_time = clock();
  // C/C++标准库读文件
  // cout<<"arrayNum"<<endl;
  // struct stat info;
  // stat(filePath, &info);
  // auto fileSize = info.st_size;
  // char *const buf = new char[fileSize];
  // FILE* fp = fopen(filePath, "r");
  // fread(buf, fileSize, 1, fp);
  // fclose(fp);
  // ifstream fin(filePath, ios::binary);
  // fin.read(&buf[0], static_cast<streamsize>(fileSize));
  // linux库读文件
  struct stat file_info;
  int fd = open(filePath, O_RDONLY);
  fstat(fd, &file_info);
  const auto fileSize = file_info.st_size;
  char *buf = (char *)mmap(NULL, fileSize, PROT_READ, MAP_PRIVATE, fd, 0);

  // const int arrayNum = fileSize / 6000;
  char *fileOffset_start_thread1 = buf;
  char *fileOffset_end_thread1 = buf + fileSize / 2;
  char *fileOffset_start_thread2 = fileOffset_end_thread1;
  char *fileOffset_end_thread2 = buf + fileSize * 2 / 3;
  char *fileOffset_start_thread3 = fileOffset_end_thread2;
  char *fileOffset_end_thread3 = buf + fileSize;
  // char *fileOffset_start_thread4 = fileOffset_end_thread3;
  // char *fileOffset_end_thread4 = buf + fileSize;
  // vector<vector<dataType>> x_input(arrayNum+1, vector<dataType>(1024));
  vector<char> predictData(fileSize / 3000, '\n');  //用于存储预测结果的buffer

  // cout<<"point2";
  thread t1(analData_thread, predictData.begin(), buf, fileOffset_end_thread1,
            weight, 0);
  thread t2(analData_thread, predictData.begin(), fileOffset_start_thread2,
            fileOffset_end_thread2, weight, 1);
  // thread t3(analData_thread, predictData.begin(), fileOffset_start_thread3,
  // fileOffset_end_thread3, weight, 2);

  analData_thread(predictData.begin(), fileOffset_start_thread3,
                  fileOffset_end_thread3, weight, 2);
  t1.join();
  t2.join();
  close(fd);
  munmap(buf, fileSize);

  ofstream fout(filePath_save, ios::binary);
  fout.write(predictData.data(), predictPos[predictNum] * 2);
}
int main() {
  predictPos[0] = 0;
  // clock_t start_time = clock();
  unique_lock<mutex> locker_weight(g_mutex_weight);
  //主线程获取权重锁
  vector<vector<dataType>> weight(3, vector<dataType>(1024));
  //权重文件3*1024，0行为w向量（有效维为1000，第1002位为bias），1行为means，2行为var，用于标准化
  // cout<<"point3";
  thread t1(readtestFile, "/data/test_data.txt", "/projects/student/result.txt",
            weight.begin());  /// projects/student/result.txt
  /*
          读测试文件，里面又分了两个子线程解析，所以共三个子线程解析，并等待weight锁，随后等待子线程在测试集位置的定位锁，最后每个子线程进行独自预测，
          等待所有子线程返回后，统一写文件
  */
  readFile("/data/train_data.txt", weight.begin());
  //主线程用于读训练数据以及训练

  isWeight = true;
  locker_weight.unlock();
  cv_weight.notify_all();
  //训练数据返回后释放权重锁
  t1.join();
}