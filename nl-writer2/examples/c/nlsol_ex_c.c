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
#include "api/c/nlsol-c.h"

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
           "  <this_exe> highs \"\" text /tmp/stub\n"
           "  <this_exe> gurobi \"mip:return_gap=1\"");
    exit(0);
  }

  const char* solver = (argc>1) ? argv[1] : "highs";
  const char* sopts = (argc>2) ? argv[2] : "";
  int binary = (argc<=3) || strcmp("text", argv[3]);
  const char* stub = (argc>4) ? argv[4] : "stub";

  // Create custom interface
  CAPIExample example = MakeCAPIExample_Linear_01();
  NLFeeder2_C feeder = MakeNLFeeder2_C(&example, binary);
  SOLHandler2_C solhnd = MakeSOLHandler2_C(&example);
  NLUtils_C utils = MakeNLUtils_C();

  // Create NLSOL_C
  NLSOL_C nlsol = NLW2_MakeNLSOL_C(&feeder, &solhnd, &utils);

  // Solve
  NLW2_SetSolver_C(&nlsol, solver);
  NLW2_SetSolverOptions_C(&nlsol, sopts);
  if (0==NLW2_Solve_C(&nlsol, stub)) {
    printf("%s\n", NLW2_GetErrorMessage_C(&nlsol));
    exit(EXIT_FAILURE);
  }
  PrintSolution_C(&example, stub);

  // Destroy API-owned objects
  NLW2_DestroyNLSOL_C(&nlsol);

  // Destroy our custom interface and example data
  DestroyNLUtils_C(&utils);
  DestroySOLHandler2_C(&solhnd);
  DestroyNLFeeder2_C(&feeder);
  DestroyCAPIExample_Linear_01(&example);

  return 0;
}
