import math
import datetime
import sys
import numpy as np


class LR:
    def __init__(self, train_file_name, test_file_name, predict_result_file_name):
        self.train_file = train_file_name
        self.predict_file = test_file_name
        self.predict_result_file = predict_result_file_name
        self.max_iters = 150
        self.rate = 0.1
        self.feats = []
        self.labels = []
        self.feats_test = []
        self.labels_predict = []
        self.param_num = 0
        self.weight = []
        self.cl = 1000

    def loadDataSet(self, file_name, label_existed_flag):
        feats = []
        labels = []
        with open(file_name, 'r') as f:
            for line in f:
                allInfo = line.strip().split(',')
                if label_existed_flag == 1:
                    labels.append(allInfo[-1])
                    feats.append(allInfo[:self.cl])
                else:
                    feats.append(allInfo[:self.cl])
        feats = np.array(feats, dtype=np.float64)
        for idx, row in enumerate(feats):
            def foo(x):
                x *= 1000
                if x < 0:
                    ans = -x
                    ans = (ans - ans % 10) / 1000
                    return -ans
                else:
                    ans = x
                    ans = (ans - ans % 10) / 1000
                    return ans
            feats[idx] = [foo(it) for it in row]
        labels = np.array(labels, dtype=np.float64)
        return feats, labels

    def savePredictResult(self):
        f = open(self.predict_result_file, 'w')
        for it in self.answer:
            f.write(str(it)+"\n")
        f.close()

    def printInfo(self):
        print(self.train_file)
        print(self.predict_file)
        print(self.predict_result_file)
        print(self.feats)
        print(self.labels)
        print(self.feats_test)
        print(self.labels_predict)

    def initParams(self):
        self.weight = np.ones((self.param_num,), dtype=np.float)

    def compute(self, recNum, param_num, feats, w):
        return self.sigmod(np.dot(feats, w))

    def error_rate(self, recNum, label, preval):
        return np.power(label - preval, 2).sum()

    def foo(self, x):
        x *= 1000
        if x < 0:
            ans = -x
            ans %= 1000
            ans = (ans - ans % 10)
            return int(ans)
        else:
            ans = x
            ans %= 1000
            ans = (ans - ans % 10)
            return int(ans)

    def Train(self, samples):
        one, zero = 0, 0
        self.Vec1 = np.zeros(1000)
        self.Vec0 = np.zeros(1000)
        with open(self.train_file, 'r') as f:
            index = 0
            for line in f:
                if index >= samples:
                    break
                index += 1
                data = line.strip().split(',')
                data = np.array(data, dtype=np.float64)
                label = int(data[-1])
                data = data[:-1]
                data = [self.foo(it) for it in data]

                if label == 1:
                    self.Vec1 += data
                    one += 1
                else:
                    self.Vec0 += data
                    zero += 1
        # self.Vec0 /= zero
        # self.Vec1 /= one
        self.Sum, self.Delta = [], []
        for x1, x2 in zip(self.Vec0, self.Vec1):
            x1 = int(x1 / zero)
            x2 = int(x2 / one)
            self.Sum.append(x1+x2)
            self.Delta.append(x1-x2)

    def Predict(self):
        self.answer = []
        with open(self.predict_file, 'r') as f:
            for line in f:
                data = line.strip().split(',')
                data = np.array(data, dtype=np.float64)
                data = [self.foo(it) for it in data]

                if data[0] >= 200:
                    # print(data[0])
                    self.answer.append(1)
                    continue
                dis = 0
                for idx in range(1000):
                    dis += (data[idx]*2-self.Sum[idx])*self.Delta[idx]
                if dis < 0:
                    self.answer.append(1)
                else:
                    self.answer.append(0)
        self.savePredictResult()


def print_help_and_exit():
    print(
        "usage:python3 main.py train_data.txt test_data.txt predict.txt [debug]")
    sys.exit(-1)


def parse_args():
    debug = False
    if len(sys.argv) == 2:
        if sys.argv[1] == 'debug':
            print("test mode")
            debug = True
        else:
            print_help_and_exit()
    return debug


if __name__ == "__main__":
    debug = parse_args()
    # train_file = "../data/train_data.txt"
    # test_file = "../data/test_data.txt"
    # predict_file = "../data/result.txt"

    train_file = "/data/train_data.txt"
    test_file = "/data/test_data.txt"
    predict_file = "/projects/student/result.txt"

    lr = LR(train_file, test_file, predict_file)
    lr.Train(589)
    lr.Predict()

    if debug:
        index = 1000
        while index <= 1000:
            # lr.train(index)
            # lr.predict()
            answer_file = "../data/answer.txt"
            f_a = open(answer_file, 'r')
            f_p = open(predict_file, 'r')
            a = []
            p = []
            lines = f_a.readlines()
            for line in lines:
                a.append(int(float(line.strip())))
            f_a.close()

            lines = f_p.readlines()
            for line in lines:
                p.append(int(float(line.strip())))
            f_p.close()
            errline = 0
            for i in range(len(a)):
                if a[i] != p[i]:
                    errline += 1
            print(len(a), errline)
            accuracy = (len(a)-errline)/len(a)
            print("index: {}, accuracy: {}".format(index, accuracy))
            index += 10
