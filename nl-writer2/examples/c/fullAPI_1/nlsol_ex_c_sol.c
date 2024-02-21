#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "nlsol_ex_c_sol.h"

/// Declare, implementation in nlsol_ex_c_nl.c
NLHeader_C CAPI_ex_Header(void* pex_void);

void OnDualSolution(
    void* p_user_data, int nvals, void* p_api_data) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  for (int i=0; i<nvals; ++i)
    pex->sol_dual_[i] = NLW2_ReadSolVal(p_api_data);
}

void OnPrimalSolution(
    void* p_user_data, int nvals, void* p_api_data) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  for (int i=0; i<nvals; ++i)
    pex->sol_primal_[i] = NLW2_ReadSolVal(p_api_data);
}

void OnObjno(void* p_user_data, int on) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  pex->objno_  = on;
}

void OnSolveCode(void* p_user_data, int sc) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  pex->solve_code_ = sc;
}

void OnIntSuffix(
    void* p_user_data, NLW2_SuffixInfo_C si, void* p_api_data) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  int nitems_max[4] = {pex->n_var, pex->n_con, pex->n_obj, 1};
  int kind = si.kind_;
  int nmax = nitems_max[kind & 3];  // {vars, cons, objs, 1}
  const char* name = si.name_;
  const char* table = si.table_;
  printf("SUFFIX '%s': kind=%d\n", name, kind);
  if (table && strlen(table))
    printf("SUFFIX TABLE:\n%s\n", table);
  int i;
  int v;
  printf("%s", "NON-ZERO VALUES:");
  while (NLW2_IntSuffixNNZ(p_api_data)) {
    NLW2_ReadIntSuffixEntry(p_api_data, &i, &v);
    printf(" (%d, %d)", i, v);
    if (i<0 || i>=nmax) {
      NLW2_ReportIntSuffixError(
            p_api_data, "bad suffix element index");
      return;
    }
  }
  if (NLW2_IntSuffixReadOK(p_api_data))    // Can check
    printf(" ...%s\n", "OK");
  else
    printf(" ...%s\n", "FAILURE");
}

void OnDblSuffix(
    void* p_user_data, NLW2_SuffixInfo_C si, void* p_api_data) {
  CAPIExample* pex = (CAPIExample*)p_user_data;
  int nitems_max[4] = {pex->n_var, pex->n_con, pex->n_obj, 1};
  int kind = si.kind_;
  int nmax = nitems_max[kind & 3];  // {vars, cons, objs, 1}
  const char* name = si.name_;
  const char* table = si.table_;
  printf("SUFFIX '%s': kind=%d\n", name, kind);
  if (table && strlen(table))
    printf("SUFFIX TABLE:\n%s\n", table);
  int i;
  double v;
  printf("%s", "NON-ZERO VALUES:");
  while (NLW2_DblSuffixNNZ(p_api_data)) {
    NLW2_ReadDblSuffixEntry(p_api_data, &i, &v);
    printf(" (%d, %g)", i, v);
    if (i<0 || i>=nmax) {
      NLW2_ReportDblSuffixError(
            p_api_data, "bad suffix element index");
      return;
    }
  }
  if (NLW2_DblSuffixReadOK(p_api_data))    // Can check
    printf(" ...%s\n", "OK");
  else
    printf(" ...%s\n", "FAILURE");
}


NLW2_SOLHandler_C MakeSOLHandler_C(CAPIExample* pex) {
  NLW2_SOLHandler_C result
      = NLW2_MakeSOLHandler_C_Default();

  result.p_user_data_ = pex;

  result.Header = CAPI_ex_Header;
  // solve message: use default (print it)
  // AMPL options: use default
  result.OnDualSolution = OnDualSolution;
  result.OnPrimalSolution = OnPrimalSolution;
  result.OnObjno = OnObjno;
  result.OnSolveCode = OnSolveCode;
  result.OnIntSuffix = OnIntSuffix;
  result.OnDblSuffix = OnDblSuffix;

  return result;
}

void DestroySOLHandler_C(NLW2_SOLHandler_C* p) {
  p->p_user_data_ = NULL;
  NLW2_DestroySOLHandler_C_Default(p);
}
