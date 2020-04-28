#include <iostream>
#define u32 uint32_t
using namespace std;

class XJBG {
 public:
  void Simulation();
  u32 HashID(const u32 &x);

 public:
  const int N = 10;
};
void XJBG::Simulation() {}
u32 XJBG::HashID(const u32 &x) {}

int main() {
  XJBG *xjbg = new XJBG();
  xjbg->Simulation();
  return 0;
}