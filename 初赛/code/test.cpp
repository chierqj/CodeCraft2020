#include <fstream>
#include <iostream>
#include <vector>

int main() {
  std::vector<std::vector<int>> ans;
  int n = 100000000;
  for (int i = 0; i < n; ++i) {
    if (i % 30 == 0) ans.emplace_back(std::vector<int>(7, -1));
  }
  std::cerr << ans.size() << "\n";
  return 0;
}