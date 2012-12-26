#include <iostream>
#include <memory>
#include "ilogcp.h"
#include "solvers/util/expr.h"
#include "asl.h"

using std::cerr;
using std::endl;

int main(int, char **argv) {
  // Solver should be destroyed after any IloException is handled.
  std::auto_ptr<ampl::IlogCPSolver> s;
  try {
    s.reset(new ampl::IlogCPSolver());
    return s->Run(argv);
  } catch (const IloException &e) {
    cerr << "Error: " << e << endl;
  } catch (const ampl::Error &e) {
    cerr << "Error: " << e.what() << endl;
  }
  return 1;
}
