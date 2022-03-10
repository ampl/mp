
#include <stdio.h>

#include "gurobi-ampls.h"
#include "gurobidirect/gurobi-ampls-c-api.h"

/// Print warnings and/or errors
void AMPLSPrintMessages(AMPLS_MP_Solver* slv);

int RunGurobiAMPLS(const char* nl_filename, const char* slv_opt) {
  AMPLS_MP_Solver slv;
  int ret=-1;
  if (!(ret = AMPLSOpenGurobi(&slv, slv_opt)))
    if (!(ret = AMPLSLoadNLModel(&slv, nl_filename))) {
      GRBmodel* mdl = GetGRBmodel(&slv);
      if (!(ret = GRBoptimize(mdl))) {  // Optimize: doing ourselves
        ret = AMPLSReportResults(&slv);
      } else {
        AMPLSAddMessage(&slv, GRBgeterrormsg(GRBgetenv(mdl)));
      }
    }

  AMPLSPrintMessages(&slv);
  AMPLSCloseGurobi(&slv);
  return ret;
}

void AMPLSPrintMessages(AMPLS_MP_Solver* slv) {
  for (const char* const* msg=AMPLSGetMessages(slv);
       NULL!=*msg; ++msg)
    printf("%s\n", *msg);
}
