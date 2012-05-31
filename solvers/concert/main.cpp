#include <iostream>
#include "concert.h"
#include "util.h"
#include "solvers/asl.h"

using std::cerr;
using std::endl;

int main(int argc, char **argv) {
  try {
    return Driver().run(argc, argv);
  } catch (const IloException &e) {
    cerr << "Error: " << e << endl;
  } catch (const Error &e) {
    cerr << "Error: " << e.what() << endl;
  }
}
