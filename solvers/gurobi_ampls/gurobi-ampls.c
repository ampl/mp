
#include <stdio.h>

#include "gurobi-ampls.h"
#include "gurobi/gurobi-ampls-c-api.h"

/// Print warnings and/or errors
void AMPLSPrintMessages(AMPLS_MP_Solver* slv);

int RunGurobiAMPLS(const char* nl_filename, const char* slv_opt) {
  AMPLS_MP_Solver* pslv;
  int ret=-1;
  if (!(pslv = AMPLSOpenGurobi(slv_opt)))
    if (!(ret = AMPLSLoadNLModel(pslv, nl_filename))) {
      GRBmodel* mdl = GetGRBmodel(pslv);
      if (!(ret = GRBoptimize(mdl))) {  // Optimize: doing ourselves
        ret = AMPLSReportResults(pslv);
      } else {
        AMPLSAddMessage(pslv, GRBgeterrormsg(GRBgetenv(mdl)));
      }
    }

  AMPLSPrintMessages(pslv);
  AMPLSCloseGurobi(pslv);
  return ret;
}

void AMPLSPrintMessages(AMPLS_MP_Solver* slv) {
  for (const char* const* msg=AMPLSGetMessages(slv);
       NULL!=*msg; ++msg)
    printf("%s\n", *msg);
}
