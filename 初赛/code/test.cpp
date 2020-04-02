#include <fstream>
#include <iostream>

int main() {
  std::ofstream fout("log.txt");
  for (int i = 0; i < 1000000; ++i) {
    fout << i << "," << i + 1 << "," << i + 2 << "," << i + 3 << "," << i + 4
         << "," << i + 5 << "," << i + 6 << "\n";
  }
  fout.close();
  return 0;
}