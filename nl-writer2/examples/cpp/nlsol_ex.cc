/**
 * Example using NLWriter2 and SOLReader2
 * to write NL file, execute solver, and read SOL file.
 *
 */

#include <cstdlib>

#include "nlsol_ex_mdl.h"
#include "nlsol_ex_nl.h"
#include "nlsol_ex_sol.h"

#include "mp/nlsol.h"
#include "mp/nl-writer2.h"
#include "mp/nl-writer2.hpp"
#include "mp/sol-reader2.h"
#include "mp/sol-reader2.hpp"

/// Invoke:
///   (this_exe) ipopt ["outlev=1 timelim=500" [text [filestub]]]
int main(int argc, const char* const* argv) {
  if (argc<2) {
    printf("%s\n",
           "AMPL NL writer example.\n"
           "Usage:\n"
           "  <this_exe> <solver> [\"<solver_options>\" [binary/text [<stub>]]],\n\n"
           "where <solver> is ipopt, gurobi, minos, ...;\n"
           "binary/text is the NL format (default: binary.)\n"
           "Examples:\n"
           "  <this_exe> ipopt \"\" text /tmp/stub\n"
           "  <this_exe> gurobi \"nonconvex=2 funcpieces=-2 funcpieceratio=1e-4\"");
    exit(0);
  }
  std::string solver = (argc>1) ? argv[1] : "minos";
  std::string sopts = (argc>2) ? argv[2] : "";
  bool binary = (argc<=3) || std::strcmp("text", argv[3]);
  std::string stub = (argc>4) ? argv[4] : "";

  ExampleModel emdl;
  ExampleNLFeeder2 nlf(emdl, binary);
  ExampleSOLHandler2 esolh(emdl);
  mp::NLUtils utils;

  mp::NLSOL nlsol(&utils);
  nlsol.SetFileStub(stub);
  if (!nlsol.LoadModel(nlf)
      || !nlsol.Solve(solver, sopts)
      || !nlsol.ReadSolution(esolh)) {
    printf("%s\n", nlsol.GetErrorMessage());
    return EXIT_FAILURE;
  } else {
    esolh.PrintSolution(stub);
  }

  return EXIT_SUCCESS;
}
