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

#if MSK_VERSION_MAJOR >= 12
namespace {

/// The values of the AMPL options for the optimizer parameters and for
/// folding follow the numbering of the MOSEK 11 enumerations, which
/// scripts depend on. MOSEK 12 renumbered the optimizer enumeration and
/// dropped a folding mode, so the values are translated.
struct ValueMap { int ampl; int mosek; };

const ValueMap optimizer_to_mosek[] = {
  {0, MSK_OPTIMIZER_CONIC},
  {1, MSK_OPTIMIZER_DUAL_SIMPLEX},
  {2, MSK_OPTIMIZER_FREE},
  {3, MSK_OPTIMIZER_FREE_SIMPLEX},
  {4, MSK_OPTIMIZER_INTPNT},
  {5, MSK_OPTIMIZER_MIXED_INT},
  {6, MSK_OPTIMIZER_DUAL_SIMPLEX},        // MOSEK 11: new dual simplex
  {7, MSK_OPTIMIZER_PRIMAL_SIMPLEX},      // MOSEK 11: new primal simplex
  {8, MSK_OPTIMIZER_PRIMAL_SIMPLEX},
  {9, MSK_OPTIMIZER_PRIMAL_DUAL_SIMPLEX},
  {10, MSK_OPTIMIZER_CONCURRENT}
};
/// Reverse direction: the deprecated interior-point optimizer, which
/// has no AMPL value, is shown as the interior-point optimizer.
const ValueMap optimizer_from_mosek[] = {
  {10, MSK_OPTIMIZER_CONCURRENT},
  {0, MSK_OPTIMIZER_CONIC},
  {4, MSK_OPTIMIZER_DEPRECATED_INTPNT},
  {1, MSK_OPTIMIZER_DUAL_SIMPLEX},
  {2, MSK_OPTIMIZER_FREE},
  {3, MSK_OPTIMIZER_FREE_SIMPLEX},
  {4, MSK_OPTIMIZER_INTPNT},
  {5, MSK_OPTIMIZER_MIXED_INT},
  {9, MSK_OPTIMIZER_PRIMAL_DUAL_SIMPLEX},
  {8, MSK_OPTIMIZER_PRIMAL_SIMPLEX}
};
const ValueMap folding_to_mosek[] = {
  {0, MSK_FOLDING_MODE_OFF},
  {1, MSK_FOLDING_MODE_FREE},
  {2, MSK_FOLDING_MODE_FREE},             // MOSEK 11: free unless basic
  {3, MSK_FOLDING_MODE_FORCE}
};
const ValueMap folding_from_mosek[] = {
  {0, MSK_FOLDING_MODE_OFF},
  {1, MSK_FOLDING_MODE_FREE},
  {3, MSK_FOLDING_MODE_FORCE}
};

template <size_t N>
int Translate(const ValueMap (&map)[N], int value, bool to_mosek) {
  for (const auto& m : map)
    if (value == (to_mosek ? m.ampl : m.mosek))
      return to_mosek ? m.mosek : m.ampl;
  return value;   // not an AMPL value: MOSEK rejects it, or it is displayed as is
}

/// Optimizer parameters share one numbering
bool IsOptimizerParam(MSKiparame key) {
  return MSK_IPAR_OPTIMIZER == key || MSK_IPAR_MIO_ROOT_OPTIMIZER == key
      || MSK_IPAR_MIO_NODE_OPTIMIZER == key;
}

int AmplToMosekValue(MSKiparame key, int value) {
  if (IsOptimizerParam(key))
    return Translate(optimizer_to_mosek, value, true);
  if (MSK_IPAR_FOLDING_USE == key)
    return Translate(folding_to_mosek, value, true);
  return value;
}

int MosekToAmplValue(MSKiparame key, int value) {
  if (IsOptimizerParam(key))
    return Translate(optimizer_from_mosek, value, false);
  if (MSK_IPAR_FOLDING_USE == key)
    return Translate(folding_from_mosek, value, false);
  return value;
}

}  // namespace
#else
namespace {
// Same numbering as the AMPL options
int AmplToMosekValue(MSKiparame, int value) { return value; }
int MosekToAmplValue(MSKiparame, int value) { return value; }
}  // namespace
#endif

void MosekCommon::GetSolverOption(MSKiparame key, int &value) const {
  MOSEK_CCALL(MSK_getintparam(lp(), key, &value));
  value = MosekToAmplValue(key, value);
}

void MosekCommon::SetSolverOption(MSKiparame key, int value) {
  // Rejection messages show the value as given by the user
  CheckParamSet(MSK_putintparam(lp(), key, AmplToMosekValue(key, value)),
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
