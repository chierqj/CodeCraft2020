#include <iostream>

int main() {
  for (int i = 0; i < 1000; ++i) {
    std::cout << "num = (*(ptr  + move + " << i * 6 + 2
              << ") - '0') * 200 + (*(ptr + move + " << i * 6 + 3
              << ") - '0') * 20;\n";
    std::cout << "sum = m_PredictSum[" << i << "];\n";
    std::cout << "delta = m_PredictDelta[" << i << "];\n";
    std::cout << "distance += (num - sum) * delta;\n";
  }
  return 0;
}