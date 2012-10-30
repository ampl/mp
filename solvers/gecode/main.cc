#include "gecode.h"

int main(int, char **argv) {
  // TODO: catch exceptions
  ampl::GecodeDriver d;
  d.run(argv);
  return 0;
}
