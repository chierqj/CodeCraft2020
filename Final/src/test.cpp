#include <unistd.h>

#include <iomanip>
#include <iostream>
#include <string>
/*设置必备的头文件*/
using namespace std;

void print_process(string strName, float fValue, bool flag = false) {
  if (0 > fValue || fValue > 1) {
    cout << "进度参数必须大于0且小于1" << endl;
    return;
  }

  fValue = fValue * 100;

  cout << setiosflags(ios::fixed) << setprecision(2);

  string tag =
      "[" + strName + "]" + "process:" + string((fValue / 10), '*') + "[";
  // flush擦除，\r定位到行首
  cout << flush << '\r' << tag << fValue << "%]";
  if (flag) {
    usleep(1000);  // 1000us
  }
}

int main() {
  for (int i = 0; i <= 12345; i++) {
    // flush擦除，\r定位到行首
    float f = (float)i / 12345.0;
    print_process("mytest", f, true);
  }
  cout << endl;
  cin.get();
  return 0;
}