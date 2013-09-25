#include "solvers/ilogcp/cplex.h"
#include "solvers/ilogcp/concert.h"

int main(int, char **argv) {
  return ampl::RunSolver<ampl::CPLEXSolver>(argv);
}
