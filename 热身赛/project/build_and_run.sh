
rm -fr test
g++ -std=c++11 -O3 main.cpp -o test -lpthread -D LOCAL_TRAIN

./test