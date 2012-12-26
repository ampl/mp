#include <iostream>
#include "solvers/gecode/gecode.h"

int main(int, char **argv) {
  try {
    return ampl::GecodeSolver().Run(argv);
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
  return 1;
}
