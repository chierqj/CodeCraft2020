#/ bin / bash
rm -fr main
g++ -std=c++11 -O3 chenchen.cpp -o chen -lpthread
time ./chen