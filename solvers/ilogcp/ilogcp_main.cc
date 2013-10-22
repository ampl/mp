#include "solvers/ilogcp/ilogcp.h"

int main(int, char **argv) {
  // Solver should be destroyed after any IloException is handled.
  std::auto_ptr<ampl::IlogCPSolver> s;
  try {
    s.reset(new ampl::IlogCPSolver());
    return s->Run(argv);
  } catch (const IloException &e) {
    std::cerr << "Error: " << e << std::endl;
  } catch (const ampl::Error &e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
  return 1;
}
