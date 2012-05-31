#include "concert.h"
#include "util.h"
#include "solvers/asl.h"

int main(int argc, char **argv) {
  try {
    return Driver().run(argc, argv);
  } catch (const Error& e) {
    Printf("%s\n", e.what());
  }
}
