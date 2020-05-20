#!/bin/bash

RED='\e[1;31m'  # 红
GREEN='\e[1;32m'  # 绿
YELLOW='\033[1;33m'  # 黄
BLUE='\E[1;34m'  # 蓝
PINK='\E[1;35m'  # 粉红
RES='\033[0m'  # 清除颜色


test_file[0]="../data/data0/test_data.txt"
test_file[1]="../data/data1/test_data.txt"
test_file[2]="../data/data2/test_data.txt"
test_file[3]="../data/data3/test_data.txt"
test_file[4]="../data/data4/test_data.txt"
test_file[5]="../data/data5/test_data.txt"
test_file[6]="../data/data6/test_data.txt"
test_file[7]="../data/data7/test_data.txt"
test_file[8]="../data/data8/test_data.txt"
test_file[9]="../data/data9/test_data.txt"
test_file[10]="../data/data10/test_data.txt"
test_file[11]="../data/data11/test_data.txt"
test_file[12]="../data/data12/test_data.txt"
test_file[13]="../data/std/test_data.txt"
test_file[14]="../data/big/test_data.txt"


answer_file[0]="../data/data0/answer.txt"
answer_file[1]="../data/data1/answer.txt"
answer_file[2]="../data/data2/answer.txt"
answer_file[3]="../data/data3/answer.txt"
answer_file[4]="../data/data4/answer.txt"
answer_file[5]="../data/data5/answer.txt"
answer_file[6]="../data/data6/answer.txt"
answer_file[7]="../data/data7/answer.txt"
answer_file[8]="../data/data8/answer.txt"
answer_file[9]="../data/data9/answer.txt"
answer_file[10]="../data/data10/answer.txt"
answer_file[11]="../data/data11/answer.txt"
answer_file[12]="../data/data12/answer.txt"
answer_file[13]="../data/std/answer.txt"
answer_file[14]="../data/big/answer.txt"




g++ -std=c++11 -O3 $1 -o main -lpthread -fpic

if [ ! -d "/data/" ]
then
    mkdir /data
fi

if [ ! -d "/projects/student/" ]
then
    mkdir /projects
    mkdir /projects/student
fi


mount -t tmpfs -o size=1g tmpfs /data
mount -t tmpfs -o size=2g tmpfs /projects/student

echo "----------------------------------"

date > log.txt
old=1
cnt=0
total_cnt=0

for ((i=0;i<${#test_file[@]};i++))
do
    echo -e "${PINK}[File$i: ${test_file[$i]}]${RES}"
    cp ${test_file[$i]} /data/test_data.txt
    time ./main
    if [ -f "/projects/student/result.txt" ]
    then
        diff /projects/student/result.txt "${answer_file[$i]}" >> log.txt
        len=$(sed -n '$=' log.txt)
        if [ $len -gt $old ]
        then
            echo -e "${RED}[@ Wrong Answer]${RES}"
        else
            echo -e "${GREEN}[@ Accepted]${RES}"
            let cnt+=1
        fi
        rm /projects/student/result.txt
    else
        echo -e "${RED}[@ File not generated${RES}"
    fi
        old=$len
    let total_cnt+=1
    echo "----------------------------------"
done

echo -e "${YELLOW}[$cnt/$total_cnt Accepted]${RES}" | tee -a ./log.txt