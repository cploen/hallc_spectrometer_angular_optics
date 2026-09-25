#include "annotation_layout.hh"
#include <iostream>
#include <string>

int main(int argc, char** argv) {
  using namespace hms_display;
  if (argc == 2 && std::string(argv[1]) == "--aliases") {
    for (const auto& q : quantities)
      std::cout << q.conceptTag << '\t' << q.code << '\t' << q.branch << '\n';
  } else if (argc == 2 && std::string(argv[1]) == "--layers") {
    for (const auto& layer : layers) std::cout << layer.name << '\n';
  } else {
    std::string name;
    while (std::cin >> name)
      for (const auto& layer : layers)
        if (primitiveLayer(name) == layer.layer)
          std::cout << name << '\t' << layer.name << '\n';
  }
}
