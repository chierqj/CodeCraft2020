# /bin/bash
rm -fr main
g++ -std=c++11 main.cpp -o main -lpthread -fpic -D LOCAL -D TESTSPEED
time ./main