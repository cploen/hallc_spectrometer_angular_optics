#include "spectrometer_config.h"
#include <iomanip>
#include <iostream>
int main() {
  const auto spec=hallc::profileForName("HMS");
  double a,z,x,y,rx,ry,bpm;
  std::cout << std::setprecision(17);
  while (std::cin >> a >> z >> x >> y >> rx >> ry >> bpm) {
    const auto t=hallc::targetTruth(spec,a,z,x,y,0,rx,ry,bpm);
    std::cout << t.xtar << ' ' << t.ytar << ' ' << t.xptar << ' ' << t.yptar << '\n';
  }
}
