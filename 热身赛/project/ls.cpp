#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>
using namespace std;

vector<int> Index;
//矩阵操作
struct Matrix {
    typedef vector<double> Mat1;
    typedef vector<vector<double> > Mat2;
};
void out(Matrix::Mat1 mat) {
    for (auto &x : mat) {
        cout << x << " ";
    }
    cout << "\n";
}
void out(Matrix::Mat2 mat) {
    for (auto &x : mat) {
        out(x);
    }
}
//点乘
static double Dot(const Matrix::Mat1 &mat1, const Matrix::Mat1 &mat2) {
    int n = mat1.size();
    double ans = 0;
    for (int i = 0; i < n; i++) {
        ans += mat1[i] * mat2[i];
    }
    return ans;
};
//乘法
static Matrix::Mat1 operator*(const Matrix::Mat2 &mat1,
                              const Matrix::Mat1 &mat2) {
    int n = mat1.size();
    Matrix::Mat1 mat;
    for (const auto &x : mat1) {
        mat.emplace_back(Dot(x, mat2));
    }
    return mat;
};
static Matrix::Mat2 operator*(const Matrix::Mat1 &mat1,
                              const Matrix::Mat1 &mat2) {
    int n = mat1.size(), m = mat2.size(), id = 0;
    Matrix::Mat2 mat(n);
    for (const auto &x : mat1) {
        for (const auto &y : mat2) {
            mat[id].emplace_back(x * y);
        }
        id++;
    }
    return mat;
};
//转置
static Matrix::Mat2 T(const Matrix::Mat2 &mat1) {
    int n = mat1.size(), m = mat1[0].size();
    Matrix::Mat2 mat(m, Matrix::Mat1(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            mat[j][i] = mat1[i][j];
        }
    }
    return mat;
};
static Matrix::Mat1 operator-(const Matrix::Mat1 &mat1,
                              const Matrix::Mat1 &mat2) {
    int n = mat1.size();
    Matrix::Mat1 mat;
    for (int i = 0; i < n; i++) {
        mat.emplace_back(mat1[i] - mat2[i]);
    }
    return mat;
}
static Matrix::Mat2 operator*(const Matrix::Mat2 &mat1,
                              const Matrix::Mat2 &mat2) {
    Matrix::Mat2 mat2T = T(mat2);
    int n = mat1.size(), m = mat2[0].size();
    Matrix::Mat2 mat(n, Matrix::Mat1(m));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            mat[i][j] = Dot(mat1[i], mat2T[j]);
        }
    }
    return mat;
};

class LR {
   public:
    LR(string trainF, string testF, string predictOutF, string answerF);
    void init();
    void predict();             //预测
    void judge(string &file);   //判分
    void getNumber(int start);  //获取训练、测试数据信息

   private:
    Matrix::Mat2 trainX;  //训练数据集
    Matrix::Mat1 trainY;  //训练集标签
    Matrix::Mat2 testX;   //测试数据集
    Matrix::Mat1 testY;   //预测的结果
    vector<int> answer;   //进行比对的答案
    Matrix::Mat1 pred;    //预测的sigmod结果
    Matrix::Mat1 weight;  //权重矩阵
    string trainFile;
    string testFile;
    string answerFile;
    string predictOutFile;

   private:
    double learningRate = 0.01;  //学习率
    int iterations = 160;        //训练轮数
    int feature = 0;             //特征数
    int trainNum = 1000;         //样本数
    int YLabTrain = 0;           //样本中正标签
    int NLabTrain = 0;           //样本中负标签
    int YLabTest = 0;            //测试结果中正标签
    int NLabTest = 0;            //测试结果中负标签
    double totTime = 0;          //运行时间

   private:
    void loadData(const string &file, Matrix::Mat2 &vt1, Matrix::Mat1 &vt2);
    void writeWeight(Matrix::Mat1 &vt, string &file);
    void writeData(Matrix::Mat1 &vt, string &file);
    void LoadChar(const vector<char> &c, Matrix::Mat2 &vt1, Matrix::Mat1 &vt2,
                  const bool &flag);
    void changeWeight(Matrix::Mat1 &vt);
    void NewloadData(const string &file, Matrix::Mat2 &vt1, Matrix::Mat1 &vt2);
    void setFeature(Matrix::Mat2 &mat1, Matrix::Mat2 &mat2);
    void init_weight();
    void train();
    int getLab(double &z);
    double sigmod(double &x);
};

double LR::sigmod(double &x) { return 1.0 / (1.0 + exp(-x)); }
LR::LR(string trainF, string testF, string predictOutF, string answerF) {
    trainFile = trainF;
    testFile = testF;
    predictOutFile = predictOutF;
    answerFile = answerF;
    init();
    train();
}

inline void ShuffleVector(std::vector<int> &vt) {
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    shuffle(vt.begin(), vt.end(), default_random_engine(seed));
}
void LR::getNumber(int start) {
    YLabTest = 0;
    YLabTrain = 0;
    NLabTest = 0;
    NLabTrain = 0;
    cout << "================\n";
    for (auto &x : trainY) {
        if (x == 1) {
            YLabTrain++;
        } else
            NLabTrain++;
    }
    cout << "================\n";
    for (auto &x : testY) {
        if (x > 0.5) {
            YLabTest++;
        } else
            NLabTest++;
    }
    int end = clock();
    totTime = (double)(end - start) / CLOCKS_PER_SEC;
    cout << "================\n";
    cout << "训练集样本数：" << trainNum << "\n";
    cout << "训练集特征数：" << feature << "\n";
    cout << "训练集正样本数：" << YLabTrain << "\n";
    cout << "训练集负样本数：" << NLabTrain << "\n";
    cout << "================\n";
    cout << "测试集样本数：" << testY.size() << "\n";
    cout << "测试集正样本数：" << YLabTest << "\n";
    cout << "测试集负样本数：" << NLabTest << "\n";
    cout << "运行时间：" << totTime << "s\n";
    cout << "================\n";
}

void LR::init_weight() {
    weight = Matrix::Mat1(feature, 1);
    weight.emplace_back(0);
    return;
    int limit = sqrt(1.0 / feature) * 10000;
    double p;
    for (int i = 0; i < feature; i++) {
        p = (rand() % (limit * 2) - limit) * 1.0 / 10000;
        weight.emplace_back(p);
    }
    weight.emplace_back(0);
}

void LR::setFeature(Matrix::Mat2 &mat1, Matrix::Mat2 &mat2) {
    Matrix::Mat1 mat3;
    for (auto &x : mat1) {
        mat3.clear();
        for (int i = 0; i < feature; i++) {
            mat3.emplace_back(x[Index[i]]);
        }
        mat3.emplace_back(1);
        mat2.emplace_back(mat3);
    }
}
void LR::init() {
    cout << "开始初始化\n";
    NewloadData(trainFile, trainX, trainY);
    loadData(testFile, testX, testY);
    trainNum = trainY.size();
    feature = trainX[0].size() - 1;
    init_weight();
}

//////////////////

double ToDouble(string &s) {
    int n = s.size(), i = n - 1;
    double ans = 0;
    if (s[0] == '-') {
        return 0;
    } else if (s[0] == '0') {
        while (i > 1) {
            ans = ans + (s[i--] - '0');
            ans /= 10;
        }
        return ans;
    } else {
        return 1;
    }
}
void LR::NewloadData(const string &file, Matrix::Mat2 &vt1, Matrix::Mat1 &vt2) {
    ifstream fin(file);
    std::vector<double> features;
    string line, temp;
    int label;
    while (fin) {
        getline(fin, line);
        if (file == trainFile && vt1.size() == trainNum) {
            break;
        }
        if (line.empty()) {
            continue;
        }
        istringstream iss(line);
        features.clear();
        while (getline(iss, temp, ',')) {
            features.emplace_back(ToDouble(temp));
        }
        if (file == trainFile) {
            label = features.back();
            features.pop_back();
            features.emplace_back(1);
            vt1.emplace_back(features);
            vt2.emplace_back(label);
        } else {
            features.emplace_back(1);
            vt1.emplace_back(features);
        }
    }
    fin.close();
}

////////////////////

void LR::LoadChar(const vector<char> &ch, Matrix::Mat2 &vt1, Matrix::Mat1 &vt2,
                  const bool &flag) {
    bool ne = false, dot = false, one = false;
    int r = 10;
    double sum = 0, x;
    Matrix::Mat1 features;
    for (auto &c : ch) {
        if (c == ',' || c == '\n') {
            if (ne) {
                sum = 0;
                ne = false;
            }
            if (one) {
                sum = 1;
                one = false;
            }
            dot = false;
            r = 10;
            if (c == ',') {
                features.emplace_back(sum);
            } else {
                if (flag) {
                    features.emplace_back(1);
                    vt1.emplace_back(features);
                    vt2.emplace_back(sum);
                    if (vt1.size() == trainNum) return;
                } else {
                    features.emplace_back(sum);
                    features.emplace_back(1);
                    vt1.emplace_back(features);
                }
                features.clear();
            }
            sum = 0;
        } else {
            if (ne) {
                continue;
            }
            if (one) {
                continue;
            }
            if (c == '-') {
                ne = true;
                continue;
            }
            if (c == '.') {
                dot = true;
                continue;
            }
            x = c - '0';
            if (!dot) {
                if (c == '0')
                    sum = 0;
                else
                    one = true;

            } else {
                sum += x / r;
                r *= 10;
            }
        }
    }
}

void LR::loadData(const string &file, Matrix::Mat2 &vt1, Matrix::Mat1 &vt2) {
    ifstream fin(file, ios::binary);
    bool flag = (file == trainFile);
    std::vector<char> buf(
        static_cast<unsigned int>(fin.seekg(0, ios::end).tellg()));
    fin.seekg(0, std::ios::beg)
        .read(&buf[0], static_cast<std::streamsize>(buf.size()));
    fin.close();
    LoadChar(buf, vt1, vt2, flag);
}

void LR::writeWeight(Matrix::Mat1 &vt, string &file) {
    std::string line;
    int i;
    std::ofstream fout(file);
    if (!fout.is_open()) {
        std::cout << "打开写入文件失败\n";
    }
    for (auto &x : vt) {
        fout << x << ",";
    }
    fout.close();
}

void LR::train() {
    cout << "================\n";
    cout << "------开始训练------\n";
    Matrix::Mat2 mat2 = T(trainX);
    Matrix::Mat1 mat;
    Matrix::Mat1 sig;
    for (int i = 0; i < iterations; i++) {
        // if (i % 20 == 0) cout << "第 " << i << " 轮迭代\n";
        learningRate = 5.0 / (i + 1) + 0.001;
        mat = trainX * weight;
        for (auto &x : mat) {
            x = sigmod(x);
        }
        sig = mat - trainY;
        sig = mat2 * sig;
        for (int j = 0; j <= feature; j++) {
            weight[j] -= sig[j] * learningRate;
        }
    }
    string modelweight = "modelweight.txt";
    writeWeight(weight, modelweight);
    cout << "----------------\n";
    cout << "----------------\n";
    cout << "------训练结束------\n";
    cout << "================\n";
}

void LR::writeData(Matrix::Mat1 &vt, string &file) {
    std::string line;
    int i;
    std::ofstream fout(file);
    if (!fout.is_open()) {
        std::cout << "打开写入文件失败\n";
    }
    for (auto &x : vt) {
        fout << (x > 0.5 ? 1 : 0) << "\n";
    }
    fout.close();
}
void LR::judge(string &file) {
    ifstream fin(file);
    int x, cor = 0, n = testY.size();
    while (fin) {
        fin >> x;
        answer.emplace_back(x);
    }
    fin.close();
    for (int i = 0; i < testY.size(); i++) {
        if ((testY[i] > 0.5 ? 1 : 0) == answer[i]) cor++;
    }
    cout << "================\n";
    cout << "预测准确率为: " << cor * 1.0 / n << "\n";
    cout << "================\n";
}
void LR::predict() {
    testY = testX * weight;
    for (auto &x : testY) {
        x = sigmod(x);
    }
    writeData(testY, predictOutFile);
}

int main(int argc, char *argv[]) {
    srand((unsigned)time(NULL));
    int start_time = clock();
    int type = 0;
    string trainFile = "../data/train_data.txt";
    string testFile = "../data/test_data.txt";
    string predictFile = "../data/result.txt";
    string answerFile = "../data/answer.txt";
    if (type) {
        trainFile = "/data/train_data.txt";
        testFile = "/data/test_data.txt";
        predictFile = "/projects/student/result.txt";
        answerFile = "/projects/student/answer.txt";
    }

    LR logist(trainFile, testFile, predictFile, answerFile);
    cout << "================\n";
    cout << "开始预测\n";
    logist.predict();
    cout << "预测完成\n";
    cout << "================\n";
    logist.getNumber(start_time);
    logist.judge(answerFile);
    return 0;
}