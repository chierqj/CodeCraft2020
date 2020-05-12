#!/bin/bash


test_file[0]="../data/std/test_data.txt"
test_file[1]="../data/15371869/test_data.txt"
test_file[2]="../data/18908526/test_data.txt"
test_file[3]="../data/19630345/test_data.txt"


answer_file[0]="../data/std/answer.txt"
answer_file[1]="../data/15371869/answer.txt"
answer_file[2]="../data/18908526/answer.txt"
answer_file[3]="../data/19630345/answer.txt"


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
    echo ${test_file[$i]}
    cp ${test_file[$i]} /data/test_data.txt
    time ./main
    if [ -f "/projects/student/result.txt" ]
    then
        diff /projects/student/result.txt "${answer_file[$i]}" >> log.txt
        len=$(sed -n '$=' log.txt)
        if [ $len -gt $old ]
        then
            echo "Fail!"
        else
            echo "Pass!"
            let cnt+=1
        fi
        rm /projects/student/result.txt
    else
        echo "File not generated !"
    fi
        old=$len
    let total_cnt+=1
    echo "----------------------------------"
done

echo "$cnt/$total_cnt Passed" | tee -a ./log.txt