/**
 * C API example using NLWriter2 and SOLReader2
 * to write NL file, execute solver, and read SOL file.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nlsol_ex_c_model.h"
#include "nlsol_ex_c_nl.h"
#include "nlsol_ex_c_sol.h"
#include "nlsol_ex_c_nlutils.h"
#include "api/c/nl-solver-c.h"

/// main()
int main(int argc, const char* const* argv) {
  if (argc<2) {
    printf("%s\n",
           "AMPL NL writer C API example.\n"
           "Usage:\n"
           "  <this_exe> <solver> [\"<solver_options>\" [binary/text [<stub>]]],\n\n"
           "where <solver> is ipopt, gurobi, minos, ...;\n"
           "binary/text is the NL format (default: binary.)\n"
           "Examples:\n"
           "  <this_exe> highs \"writeprob=/tmp/stub.lp\" text /tmp/stub\n"
           "  <this_exe> gurobi \"mip:return_gap=1\"");
    exit(0);
  }

  int result = EXIT_SUCCESS;
  const char* solver = (argc>1) ? argv[1] : "highs";
  const char* sopts = (argc>2) ? argv[2] : "";
  int binary = (argc<=3) || strcmp("text", argv[3]);
  const char* stub = (argc>4) ? argv[4] : "";

  // Create custom interface
  CAPIExample example = MakeCAPIExample_Linear_01();
  NLW2_NLFeeder2_C feeder = MakeNLFeeder2_C(&example, binary);
  NLW2_SOLHandler2_C solhnd = MakeSOLHandler2_C(&example);
  NLW2_NLUtils_C utils = MakeNLUtils_C();

  // Create NLSOL_C
  NLW2_NLSolver_C nlsol = NLW2_MakeNLSolver_C(&utils);

  // Solve
  NLW2_SetFileStub_C(&nlsol, stub);
  if (!NLW2_SolveFeederHandler_C(&nlsol,
                                 &feeder, &solhnd,
                                 solver, sopts)) {
    printf("%s\n", NLW2_GetErrorMessage_C(&nlsol));
    result = EXIT_FAILURE;
  } else {
    PrintSolution_C(&example, stub);
  }

  // Destroy API-owned objects
  NLW2_DestroyNLSolver_C(&nlsol);

  // Destroy our custom interface and example data
  DestroyNLUtils_C(&utils);
  DestroySOLHandler2_C(&solhnd);
  DestroyNLFeeder2_C(&feeder);
  DestroyCAPIExample_Linear_01(&example);

  return result;
}
