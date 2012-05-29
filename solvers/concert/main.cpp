#include "concert.h"
#include "util.h"
#include "solvers/asl.h"

int main(int argc, char **argv) {
  try {
    return concert_main(argc, argv);
  } catch (const Error& e) {
    Printf("%s\n", e.what());
  }
}
