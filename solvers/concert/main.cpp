#include <iostream>
#include <memory>
#include "concert.h"
#include "util.h"
#include "solvers/asl.h"

using std::cerr;
using std::endl;

int main(int argc, char **argv) {
  // Driver should be destroyed after any IloException is handled.
  std::auto_ptr<Driver> d;
  try {
    d.reset(new Driver());
    return d->run(argc, argv);
  } catch (const IloException &e) {
    cerr << "Error: " << e << endl;
  } catch (const Error &e) {
    cerr << "Error: " << e.what() << endl;
  }
  return 1;
}
