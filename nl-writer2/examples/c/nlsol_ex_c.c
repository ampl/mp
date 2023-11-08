/**
 * C API example using NLWriter2 and SOLReader2
 * to write NL file, execute solver, and read SOL file.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "api/c/nl-feeder2.h"

#include "nlsol_ex_c_model.h"

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


  return 0;
}
