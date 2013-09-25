#include "solvers/ilogcp/ilogcp.h"
#include "solvers/ilogcp/concert.h"

int main(int, char **argv) {
  return ampl::RunSolver<ampl::IlogCPSolver>(argv);
}
