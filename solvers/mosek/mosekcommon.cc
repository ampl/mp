#include <cstdio>

#include "mp/format.h"
#include "mosekcommon.h"


namespace mp {

int MosekCommon::getIntAttr(MSKiinfitem_enum name)  const {
  int value = 0;
  MOSEK_CCALL(MSK_getintinf(lp(), name, &value));
  return value;
}
long long MosekCommon::getLongAttr(MSKliinfitem_enum name)  const {
  int64_t value = 0;
  MOSEK_CCALL(MSK_getlintinf(lp(), name, &value));
  return value;
}
double MosekCommon::getDblAttr(MSKdinfitem_enum name) const  {
  double value = 0;
  MOSEK_CCALL(MSK_getdouinf(lp(), name, &value));
  return value;
}

int MosekCommon::NumLinCons() const {
	return getIntAttr(MSK_IINF_OPT_NUMCON);
}

int MosekCommon::NumVars() const {
  return getIntAttr(MSK_IINF_OPT_NUMVAR);
}

int MosekCommon::NumObjs() const {
  // TODO Get number of objectives using solver API
  //return MOSEKgetnumobjs (env_, lp_);
  return 0;
}

int MosekCommon::NumQPCons() const {
  // TODO Get number of quadratic constraints using solver API
  // return getIntAttr(MOSEK_INTATTR_QCONSTRS);
  return 0;
}

int MosekCommon::NumSOSCons() const {
  // TODO Get number of SOS constraints using solver API
  // return getIntAttr(MOSEK_INTATTR_SOSS);
  return 0;
}

int MosekCommon::NumIndicatorCons() const {
  // TODO Get number of indicator constraints using solver API
  // return getIntAttr(MOSEK_INTATTR_INDICATORS);
  return 0;
}

std::string MosekCommon::ParamName(MSKparametertypee type, int key) const {
  char name[MSK_MAX_STR_LEN];
  if (MSK_RES_OK != MSK_getparamname(lp(), type, key, name))
    std::snprintf(name, sizeof(name), "parameter %d", key);
  auto it = optNames_.find({type, key});
  return it == optNames_.end()
      ? std::string(name) : fmt::format("{} (MOSEK parameter {})", it->second, name);
}

void MosekCommon::CheckParamSet(MSKrescodee code, MSKparametertypee type,
                                int key, const std::string& value) const {
  if (MSK_RES_OK == code)
    return;
  char symb[MSK_MAX_STR_LEN], desc[MSK_MAX_STR_LEN];
  MSK_getcodedesc(code, symb, desc);
  MP_RAISE(fmt::format("Rejecting value {} for {}: {} ({}); "
                       "see the MOSEK message above or docs.mosek.com "
                       "for the accepted values",
                       value, ParamName(type, key), desc, symb));
}

void MosekCommon::GetSolverOption(MSKiparame key, int &value) const {
  MOSEK_CCALL(MSK_getintparam(lp(), key, &value));
}

void MosekCommon::SetSolverOption(MSKiparame key, int value) {
  CheckParamSet(MSK_putintparam(lp(), key, value),
                MSK_PAR_INT_TYPE, key, std::to_string(value));
}

void MosekCommon::GetSolverOption(MSKdparame key, double &value) const {
  MOSEK_CCALL(MSK_getdouparam(lp(), key, &value));
}

void MosekCommon::SetSolverOption(MSKdparame key, double value) {
  CheckParamSet(MSK_putdouparam(lp(), key, value),
                MSK_PAR_DOU_TYPE, key, fmt::format("{}", value));
}

void MosekCommon::GetSolverOption(MSKsparame key, std::string &value) const {
  char *buffer;
  auto ret= MSK_getstrparamal(lp(), key, 0, &buffer);
  if (ret == MSK_RES_OK)
    value = buffer;
  MSK_freetask(lp(), buffer);
  if (ret != MSK_RES_OK)
    MP_RAISE(fmt::format("Error {} while getting parameter {}", ret, key));

}

void MosekCommon::SetSolverOption(MSKsparame key, const std::string& value) {
  CheckParamSet(MSK_putstrparam(lp(), key, value.data()),
                MSK_PAR_STR_TYPE, key, value);
}

void handleError(MSKrescodee e, const char* fname) {
	char symb[MSK_MAX_STR_LEN];
	char str[MSK_MAX_STR_LEN];
	MSK_getcodedesc(e, symb, str);
	MSKrescodetypee et;
	MSKrescodee e2 = MSK_getresponseclass(e, &et);
	if (e2 != MSK_RES_OK) MP_RAISE(
		fmt::format(
          "Error {}({}): {}\n Error in getresponseclass: {}",
      symb, e, str, e2));
  if ((int)et >= (int)MSK_RESPONSE_TRM)
    MP_RAISE(
          fmt::format(
            "Error type {}, {}({}): {}", et, symb, e, str));
	fmt::print("\nWarning {}({}): {}\n While calling '{}'\n", symb, e, str, fname);
}

} // namespace mp
