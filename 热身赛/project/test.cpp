#include <iostream>

int main() {
  for (int i = 0; i < 64; ++i) {
    std::cout << "threadSum[pid][label][i + " << i << "] += features[i + " << i
              << "];\n";
  }
  return 0;
}