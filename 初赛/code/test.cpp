#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

int main() {
  std::vector<std::pair<int, int>> vt;
  vt.emplace_back(std::make_pair(1, 1));
  vt.emplace_back(std::make_pair(1, 2));
  vt.emplace_back(std::make_pair(1, 3));
  vt.emplace_back(std::make_pair(2, 1));
  vt.emplace_back(std::make_pair(2, 2));
  vt.emplace_back(std::make_pair(2, 3));

  std::pair<int, int> p1 = std::make_pair(1, 0);
  std::pair<int, int> p2 = std::make_pair(1, 1000);
  auto b = std::lower_bound(vt.begin(), vt.end(), p1);
  auto e = std::upper_bound(vt.begin(), vt.end(), p2);
  std::cerr << e - vt.begin() << "\n";
  for (auto &it = b; it != e; ++it) {
    std::cerr << it->first << ", " << it->second << "\n";
  }
  return 0;
}