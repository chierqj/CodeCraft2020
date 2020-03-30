import math
import datetime
import sys
import numpy as np

# train_file = "../data/train_data.txt"
# test_file = "../data/test_data.txt"
# predict_file = "../data/result.txt"

train_file = "/data/train_data.txt"
test_file = "/data/test_data.txt"
predict_file = "/projects/student/result.txt"

answer = []
Vec1 = np.zeros(1000)
Vec0 = np.zeros(1000)
Sum, Delta = [], []


def savePredictResult():
    f = open(predict_file, 'w')
    for it in answer:
        f.write(str(it)+"\n")
    f.close()


def foo(x):
    x *= 1000
    if x < 0:
        ans = -x
        ans %= 1000
        ans = (ans - ans % 10)
        return int(-ans)
    else:
        ans = x
        ans %= 1000
        ans = (ans - ans % 10)
        return int(ans)


def Train(samples):
    one, zero = 0, 0
    global Vec0, Vec1
    global Sum, Delta

    with open(train_file, 'r') as f:
        index = 0
        for line in f:
            if index >= samples:
                break
            index += 1
            data = line.strip().split(',')
            data = np.array(data, dtype=np.float64)
            label = int(data[-1])
            data = data[:-1]
            data = [foo(it) for it in data]

            if label == 1:
                Vec1 += data
                one += 1
            else:
                Vec0 += data
                zero += 1
            # if zero > 234 and zero / (one + zero) >= 0.32:
            #     break

    print(zero, one, one+zero)
    for x1, x2 in zip(Vec0, Vec1):
        x1 = int(x1 / zero)
        x2 = int(x2 / one)
        Sum.append(x1+x2)
        Delta.append(x1-x2)


def Predict():
    global Sum, Delta
    print(len(Sum), len(Delta))
    with open(test_file, 'r') as f:
        for line in f:
            data = line.strip().split(',')
            data = np.array(data, dtype=np.float64)
            data = [foo(it) for it in data]

            if data[0] >= 200:
                answer.append(1)
                continue
            dis = 0
            for idx in range(1000):
                dis += (data[idx]*2-Sum[idx])*Delta[idx]
            if dis < 0:
                answer.append(1)
            else:
                answer.append(0)
    savePredictResult()


if __name__ == "__main__":
    Train(1189)
    Predict()

    # if True:
    #     answer_file = "../data/answer.txt"
    #     f_a = open(answer_file, 'r')
    #     f_p = open(predict_file, 'r')
    #     a = []
    #     p = []
    #     lines = f_a.readlines()
    #     for line in lines:
    #         a.append(int(float(line.strip())))
    #     f_a.close()

    #     lines = f_p.readlines()
    #     for line in lines:
    #         p.append(int(float(line.strip())))
    #     f_p.close()
    #     errline = 0
    #     for i in range(len(a)):
    #         if a[i] != p[i]:
    #             errline += 1
    #     print(len(a), errline)
    #     accuracy = (len(a)-errline)/len(a)
    #     print("accuracy: {}".format(accuracy))
