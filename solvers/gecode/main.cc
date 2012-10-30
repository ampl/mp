#include <iostream>
#include "solvers/gecode/gecode.h"

int main(int, char **argv) {
  try {
    ampl::GecodeDriver d;
    return d.run(argv);
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
  return 1;
}
