#include <arm_neon.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <unistd.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <thread>
#include <vector>

//#define  TEST
#define thirdnumber 4
#define paraper 0.05
#define FeaturesNumber 1000
#define predictTrueThresh 0

using namespace std;

struct NormalParam {
  int miu[1024];
};
struct impactpara {
  int impact;
  int id;
};
struct Param {
  NormalParam NormalSetP1;
  NormalParam NormalSetP0;
  impactpara impact[1024];
  int P1;
  int P0;
  int P1_1;
  int P0_1;
};
struct Dataset {
  int features[1024];
  int label;
};
Dataset ttrainset[thirdnumber][65536];
Dataset trainDataSet[thirdnumber][8192];
int predictVec[thirdnumber][65536];

class bys {
 public:
  void predict(int *testseti, int threadnum);
  int loadModel();
  int storeModel();
  bool Trainjoin();
  bool TestDatajoin();
  bool sortimpact();
  bys(string trainFile, string testFile, string predictOutFile);

 private:
  Param param;
  string trainFile;
  string testFile;
  string predictOutFile;
  string weightParamFile = "modelweight.txt";

 private:
  bool init();
  int storePredict();
  void initParam();
  int calmn();
  bool multiload(int threadnum);
  bool TrainDatajoin();
  void multitrain(int threadnum);
  bool multiloadtest(int threadnum, char *bufferii);

 private:
  int filesize;
  int sigVal0;
  int testfilesize;
  int setnum[thirdnumber];
  int trainsetnumi[thirdnumber];
  int trainsetnum;
  const int miuInitV = 0.;
};

bys::bys(string trainF, string testF, string predictOutF) {
  trainFile = trainF;
  testFile = testF;
  predictOutFile = predictOutF;
  calmn();
  init();
}

void bys::initParam() {
  int i;
  for (i = 0; i < FeaturesNumber; i++) {
    param.impact[i].id = i;
    param.NormalSetP0.miu[i] = 0;
    param.NormalSetP1.miu[i] = 0;
  }
  int sumtemp = 0;
  for (i = 0; i < trainsetnum; i++) {
    sumtemp += trainDataSet[0][i].label;
  }
  param.P1 = sumtemp;
  param.P0 = trainsetnum - sumtemp;
  param.P1_1 = 1. / param.P1;
  param.P0_1 = 1. / param.P0;
  sigVal0 = log(param.P1) - log(param.P0);
}

bool bys::init() {
  bool status = TrainDatajoin();
  if (status != true) {
    return false;
  }
  initParam();
  return true;
}

void bys::multitrain(int threadnum) {
  int muntemp = FeaturesNumber / thirdnumber;
  int buffersum = muntemp;                 // buffer´óÐ¡
  int buffermove = threadnum * buffersum;  //Æ«ÒÆÁ¿
  if (threadnum == thirdnumber - 1)
    buffersum =
        FeaturesNumber - buffermove;  //×îºóÒ»¸öÏß³Ì¿ÉÄÜºÍÆäËûÏß³ÌÊýÄ¿²»Í¬

  int i, j;

  for (j = 0; j < trainsetnum; j++) {
    if (trainDataSet[0][j].label == 1) {
      for (i = buffermove; i < buffermove + buffersum; i++) {
        param.NormalSetP1.miu[i] +=
            trainDataSet[0][j].features[i];  //¼ÆËãP0 ºÍ P1 Ìõ¼þÏÂµÄ¾ùÖµ
      }
    } else {
      for (i = buffermove; i < buffermove + buffersum; ++i) {
        param.NormalSetP0.miu[i] +=
            trainDataSet[0][j].features[i];  //¼ÆËãP0 ºÍ P1 Ìõ¼þÏÂµÄ¾ùÖµ
      }
    }
  }
  for (i = buffermove; i < buffermove + buffersum; ++i) {
    param.NormalSetP0.miu[i] /= param.P0;  //Ìõ¼þ¸ÅÂÊP0¾ùÖµ
    param.NormalSetP1.miu[i] /= param.P1;  //Ìõ¼þ¸ÅÂÊP1¾ùÖµ
  }
}
bool bys::Trainjoin() {
  vector<std::thread> vec_threads;
  for (int i = 0; i < thirdnumber; ++i) {
    std::thread th(&bys::multitrain, this, i);
    vec_threads.emplace_back(std::move(th));  // push_back() is also OK
  }
  auto it = vec_threads.begin();
  for (; it != vec_threads.end(); ++it) {
    (*it).join();
  }
  return 1;
}

char *testbuffer;
bool bys::multiloadtest(int threadnum, char *bufferii) {
  ///////////////////////
  //¼ÆËãÆ«ÒÆÁ¿
  int muntemp = testfilesize / thirdnumber;
  int buffersum = muntemp;                 // buffer´óÐ¡
  int buffermove = threadnum * buffersum;  //Æ«ÒÆÁ¿
  if (threadnum == thirdnumber - 1)
    buffersum = testfilesize - buffermove;  //×îºóÒ»¸öÏß³Ì¿ÉÄÜºÍÆäËûÏß³ÌÊýÄ¿²»Í¬
  int setnumi = buffersum / 6000;  //
                                   //////////////////////////
  int i = 0;
  char *buffer = bufferii + buffermove;

  int features[1024];
  int j = 0;
  for (int ii = 0; ii < setnumi; ii++) {
    j = 6000 * ii;
    int sigVal = 0.;
    int predictVal = 1;
    int32x4_t sum_vec32 = vdupq_n_s32(0);
    for (i = 0; i < FeaturesNumber; i += 4) {
      int32x4_t syms0 = vdupq_n_s32('0');  //'0'
      int temp[4];
      temp[0] = buffer[j];
      temp[1] = buffer[j + 6];
      temp[2] = buffer[j + 12];
      temp[3] = buffer[j + 18];
      int32x4_t _numtemp = vld1q_s32(temp);
      int32x4_t _num = vsubq_s32(_numtemp, syms0);
      int32x4_t symsk = vdupq_n_s32(1000);
      int32x4_t sum_vec = vmulq_s32(_num, symsk);
      j += 2;
      temp[0] = buffer[j];
      temp[1] = buffer[j + 6];
      temp[2] = buffer[j + 12];
      temp[3] = buffer[j + 18];
      _numtemp = vld1q_s32(temp);
      _num = vsubq_s32(_numtemp, syms0);
      symsk = vdupq_n_s32(100);  //'0'
      sum_vec = vmlaq_s32(sum_vec, _num, symsk);
      j++;
      temp[0] = buffer[j];
      temp[1] = buffer[j + 6];
      temp[2] = buffer[j + 12];
      temp[3] = buffer[j + 18];
      _numtemp = vld1q_s32(temp);
      _num = vsubq_s32(_numtemp, syms0);
      symsk = vdupq_n_s32(10);  //'0'
      sum_vec = vmlaq_s32(sum_vec, _num, symsk);
      j++;
      temp[0] = buffer[j];
      temp[1] = buffer[j + 6];
      temp[2] = buffer[j + 12];
      temp[3] = buffer[j + 18];
      _numtemp = vld1q_s32(temp);
      _num = vsubq_s32(_numtemp, syms0);
      int32x4_t _testset = vaddq_s32(sum_vec, _num);
      int para[4];  /// param.NormalSetP0.miu[i]
      para[0] = param.NormalSetP0.miu[i];
      para[1] = param.NormalSetP0.miu[i + 1];
      para[2] = param.NormalSetP0.miu[i + 2];
      para[3] = param.NormalSetP0.miu[i + 3];
      int32x4_t _para = vld1q_s32(para);
      int32x4_t delta = vsubq_s32(
          _testset,
          _para);  //////////////(testseti[i] - param.NormalSetP0.miu[i])
      delta = vmulq_s32(delta, delta);  //
      sum_vec32 = vaddq_s32(sum_vec32, delta);
      ///////////////p1
      para[0] = param.NormalSetP1.miu[i + 0];
      para[1] = param.NormalSetP1.miu[i + 1];
      para[2] = param.NormalSetP1.miu[i + 2];
      para[3] = param.NormalSetP1.miu[i + 3];
      //_testset = vld1q_s32(temptestseti);
      _para = vld1q_s32(para);
      delta = vsubq_s32(
          _testset,
          _para);  //////////////(testseti[i] - param.NormalSetP0.miu[i])
      delta = vmulq_s32(
          delta,
          delta);  //(testseti[i] - param.NormalSetP0.miu[i])*(testseti[i] -
                   // param.NormalSetP0.miu[i])
                   ////(param.NormalSetP0.sigma[i]
      sum_vec32 = vsubq_s32(sum_vec32, delta);

      j += 20;
    }
    sigVal = vaddvq_s32(sum_vec32);
    predictVal = sigVal >= predictTrueThresh ? 1 : 0;
    predictVec[threadnum][ii] = predictVal;  //
  }
  setnum[threadnum] = setnumi;
  return 1;
}

bool bys::TestDatajoin() {
  int fd = open(testFile.c_str(), O_RDONLY);
  char *buffertest =
      (char *)mmap(NULL, testfilesize, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);
  vector<std::thread> vec_threads;
  for (int i = 0; i < thirdnumber; ++i) {
    std::thread th(&bys::multiloadtest, this, i, buffertest);
    vec_threads.emplace_back(std::move(th));  // push_back() is also OK
  }
  auto it = vec_threads.begin();
  for (; it != vec_threads.end(); ++it) {
    (*it).join();
  }
  storePredict();
  return 1;
}
bool loadAnswerData(string awFile, vector<int> &awVec) {
  ifstream infile(awFile.c_str());
  if (!infile) {
    cout << "´ò¿ª´ð°¸ÎÄ¼þÊ§°Ü" << endl;
    exit(0);
  }

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

int bys::storePredict() {
  FILE *pFile;
  int i, j;
  pFile = fopen(predictOutFile.c_str(), "w");
  if (pFile == NULL) {
    return 0;
  }
  for (i = 0; i < thirdnumber; i++) {
    for (j = 0; j < setnum[i]; j++) {
      char temp = predictVec[i][j] + '0';
      fwrite(&temp, 1, 1, pFile);
      fwrite("\n", 1, 1, pFile);
    }
  }
  fclose(pFile);
  return 0;
}

bool bys::multiload(int threadnum) {
  ///////////////////////
  //¼ÆËãÆ«ÒÆÁ¿
  int muntemp = filesize / thirdnumber;
  int buffersum = muntemp;                 // buffer´óÐ¡
  int buffermove = threadnum * buffersum;  //Æ«ÒÆÁ¿
  if (threadnum == thirdnumber - 1)
    buffersum = filesize - buffermove;  //×îºóÒ»¸öÏß³Ì¿ÉÄÜºÍÆäËûÏß³ÌÊýÄ¿²»Í¬
                                        ////////////////////////////
  int i = 0;

  FILE *pFile;
  size_t len;
  ssize_t read;
  char *line = NULL;
  pFile = fopen(trainFile.c_str(), "rb");
  if (pFile == NULL) {
    return 0;
  }
  fseek(pFile, buffermove, SEEK_SET);
  int nf = 0, ns = 0;
  unsigned int sumreal = getline(&line, &len, pFile);

  while (((read = getline(&line, &len, pFile)) != -1)) {
    sumreal = sumreal + read;
    if (sumreal > buffersum) {
      break;
    }
    int tempa = line[read - 2] - '0';
    trainDataSet[threadnum][ns].label = tempa;
    nf = 0;

    int j = 0;
    for (i = 0; i < FeaturesNumber; i += 4) {
      int32x4_t syms0 = vdupq_n_s32('0');  //'0'
      int s[4] = {1, 1, 1, 1};
      int par0[4], par1[4], par2[4], par3[4];
      for (int ii = 0; ii < 4; ++ii) {
        if (line[j] == '-') {
          s[ii] = -1;
          ++j;
        }
        par0[ii] = line[j];
        j += 2;
        par1[ii] = line[j];
        ++j;
        par2[ii] = line[j];
        ++j;
        par3[ii] = line[j];
        j += 2;
      }
      // 1
      int32x4_t _numtemp = vld1q_s32(par0);
      int32x4_t _num = vsubq_s32(_numtemp, syms0);
      int32x4_t symsk = vdupq_n_s32(1000);
      int32x4_t sum_vec = vmulq_s32(_num, symsk);
      // 2
      _numtemp = vld1q_s32(par1);
      _num = vsubq_s32(_numtemp, syms0);
      symsk = vdupq_n_s32(100);  //'0'
      sum_vec = vmlaq_s32(sum_vec, _num, symsk);
      // 3
      _numtemp = vld1q_s32(par2);
      _num = vsubq_s32(_numtemp, syms0);
      symsk = vdupq_n_s32(10);  //'0'
      sum_vec = vmlaq_s32(sum_vec, _num, symsk);
      // 4
      _numtemp = vld1q_s32(par3);
      _num = vsubq_s32(_numtemp, syms0);
      sum_vec = vaddq_s32(sum_vec, _num);

      int32x4_t symski = vld1q_s32(s);
      sum_vec = vmulq_s32(sum_vec, symski);
      vst1q_s32(trainDataSet[threadnum][ns].features + i, sum_vec);
    }
    ns++;
  }
  fclose(pFile);
  trainsetnumi[threadnum] = ns;
  return 1;
}
bool bys::TrainDatajoin() {
  vector<std::thread> vec_threads;
  for (int i = 0; i < thirdnumber; ++i) {
    std::thread th(&bys::multiload, this, i);
    vec_threads.emplace_back(std::move(th));  // push_back() is also OK
  }
  auto it = vec_threads.begin();
  for (; it != vec_threads.end(); ++it) {
    (*it).join();
  }

  for (int i = 1; i < thirdnumber; ++i) {
    memcpy(trainDataSet[0] + trainsetnumi[0], trainDataSet[i],
           (trainsetnumi[i]) * sizeof(Dataset));
    trainsetnumi[0] = trainsetnumi[0] + trainsetnumi[i];
  }

  trainsetnum = trainsetnumi[0];
  return 1;
}
int bys::calmn() {
  ifstream infile(trainFile.c_str());
  string lineTitle;
  if (!infile) {
    cout << "´ò¿ª²âÊÔÎÄ¼þÊ§°Ü" << endl;
    exit(0);
  }
  infile.seekg(0, ios::end);
  filesize = infile.tellg();
  filesize = paraper * filesize;
  infile.close();

  ifstream intestfile(testFile.c_str());
  if (!intestfile) {
    cout << "´ò¿ª²âÊÔÎÄ¼þÊ§°Ü" << endl;
    exit(0);
  }
  intestfile.seekg(0, ios::end);
  testfilesize = intestfile.tellg();
  intestfile.close();
  return true;
}
int main() {
  vector<int> answerVec;
  vector<int> predictVec;
  int correctCount;
  float accurate;
  string trainFile = "/data/train_data.txt";
  string testFile = "/data/test_data.txt";
  string predictFile = "/projects/student/result.txt";
  string answerFile = "/projects/student/answer.txt";

  bys NaiveBayes(trainFile, testFile, predictFile);
  NaiveBayes.Trainjoin();

#ifdef TEST
  cout << "ready to load answer data" << endl;
  loadAnswerData(answerFile, answerVec);
#endif

  // NaiveBayes.sortimpact();

  NaiveBayes.TestDatajoin();

#ifdef TEST
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

  accurate = ((float)correctCount) / answerVec.size();
  cout << "the prediction accuracy is " << accurate << endl;
#endif

  return 0;
}
