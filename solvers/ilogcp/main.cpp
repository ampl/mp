#include <iostream>
#include <memory>
#include "solvers/util/util.h"
#include "ilogcp.h"
#include "asl.h"

using std::cerr;
using std::endl;

int main(int, char **argv) {
  // Driver should be destroyed after any IloException is handled.
  std::auto_ptr<ampl::Driver> d;
  try {
    d.reset(new ampl::Driver());
    return d->run(argv);
  } catch (const IloException &e) {
    cerr << "Error: " << e << endl;
  } catch (const ampl::Error &e) {
    cerr << "Error: " << e.what() << endl;
  }
  return 1;
}
