#/ bin / bash
rm -fr main
g++ -std=c++11 -O3 main.cpp -o main -lpthread -D LOCAL
time ./main