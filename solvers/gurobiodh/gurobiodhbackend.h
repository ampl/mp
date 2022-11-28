#ifndef MP_GUROBIODH_BACKEND_H_
#define MP_GUROBIODH_BACKEND_H_

#include <string>
#include "gurobiodhcommon.h"
#include "gurobi/gurobibackend.h"


namespace mp {

  class GurobiODHBackend : public GurobiBackend, GurobiODHCommonInfo {
  
    /// Chance for the Backend to init solver environment, etc
    void InitOptionParsing() override;

    void OpenODH();

    void Solve() override;
  
  };
}
#endif