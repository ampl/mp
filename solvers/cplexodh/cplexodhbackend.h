#ifndef MP_GUROBIODH_BACKEND_H_
#define MP_GUROBIODH_BACKEND_H_

#include <string>
#include "odh/odhcommon.h"
#include "cplexmp/cplexmpbackend.h"


namespace mp {

  class CplexODHBackend : public CplexBackend, ODHCommonInfo {
  
    void InitOptionParsing() override;

    void Solve() override;

    ~CplexODHBackend();

    void ReportResults() override; 

    void ReportODHResults();
  };
}
#endif