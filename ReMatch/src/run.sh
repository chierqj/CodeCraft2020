# /bin/bash
rm -fr main
g++ -std=c++17 -O3 main.cpp -o main -lpthread -fpic -D LOCAL -D TESTSPEED -D LOCAL_TEST
time ./main