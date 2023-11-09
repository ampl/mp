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
#include "api/c/nl-writer2-misc-c.h"
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

  CAPIExample example = MakeCAPIExample_Linear_01();
  NLFeeder2_C feeder = MakeNLFeeder2_C();
  SOLHandler2_C handler = MakeSOLHandler2_C();
  NLUtils_C utils = NLW2_MakeNLUtils_C_Default();

  NLSOL_C nlsol = NLW2_MakeNLSOL_C(&feeder, &handler, &utils);

  // Destroy API-owned objects
  NLW2_DestroyNLSOL_C(&nlsol);

  // Destroy our example data
  DestroyCAPIExample_Linear_01(&example);

  return 0;
}
