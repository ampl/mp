#ifndef MP_GUROBIODH_BACKEND_H_
#define MP_GUROBIODH_BACKEND_H_

#include <string>
#include "odh/odhcommon.h"
#include "gurobi/gurobibackend.h"


namespace mp {

  class GurobiODHBackend : public GurobiBackend, ODHCommonInfo {
  
    void InitOptionParsing() override;

    void Solve() override;

    ~GurobiODHBackend();


    void ReportResults() override;

    void ReportGurobiResults();

  };
}
#endif