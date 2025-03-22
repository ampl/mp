#include "mp/format.h"
#include "highsmpcommon.h"

namespace mp {

  void AccObjectives::setInHighs(void* highs) const {
    if (senses.size() == 0) return;
    // Only to be used when adding a quadratic objective
    assert(senses.size()==1);
    HIGHS_CCALL(parent_->Highs_changeColsCostByRange(highs, 0, coeffs.size()-1, coeffs.data()));
    HIGHS_CCALL(parent_->Highs_changeObjectiveSense(highs, senses[0]));
  }
  void AccObjectives::setAllInHighs(void* highs) const {
    if(senses.size() == 0)
      return; // no objectives yet, like when only quad objs
    if(senses.size()==1)
      setInHighs(highs);
    else {
      int nobj = senses.size();

      std::vector<double> zeroes(nobj, 0.0);
      std::vector<double> ones(nobj, 1.0);
      std::vector<double> s(nobj, 1e-5);
      std::vector<int> p(nobj, 1);
      const double* abs = abstol.empty()  ? s.data() : abstol.data();
      const double* rel = reltol.empty()  ? s.data() : reltol.data();
      const int* pri = priority.empty()   ? p.data() : priority.data();
      const double* w = weight.empty()    ? ones.data() : weight.data();

      parent_->Highs_passLinearObjectives(highs, senses.size(),
        w, zeroes.data(), coeffs.data(), 
        abs, rel, pri);
    }
  }

  void AccObjectives::setWeights(ArrayRef<double> w) {
    weight.insert(weight.begin(), w.begin(), w.end());
  }
  void AccObjectives::setOffsets(ArrayRef<double> o) {
    offset.insert(offset.begin(), o.begin(), o.end());
  }
  void AccObjectives::setRelTols(ArrayRef<double> r) {
    reltol.insert(reltol.begin(), r.begin(), r.end());
  }
  void AccObjectives::setAbsTols(ArrayRef<double> r) {
    abstol.insert(abstol.begin(), r.begin(), r.end());
  }
  void AccObjectives::setPriorities(ArrayRef<int> p) {
    priority.insert(priority.begin(), p.begin(), p.end());
  }

  void HighsCommon::LoadHighsLibrary(bool gpu) {
    // Create library loader
    setLoader(std::make_shared<HighsLoader>());
    bool libLoaded = loader().load(mp::HighsLoader::getHighsLibraryName(gpu));
    if (!libLoaded)
      throw std::runtime_error(fmt::format("Problems loading HiGHS library:\n{}", mp::HighsLoader::getHighsLibraryName(gpu)));
}
void HighsCommon::OpenSolver() {
  int status = 0;
  void* prob = loader().Highs_create();
  set_lp(prob); 
  // Create objective accumulator after the loader, as it
  // needs the reference to it
  setObjectiveAccumulator(std::make_shared<AccObjectives>());
  if (status)
    throw std::runtime_error( fmt::format(
          "Failed to create problem, error code {}.", status ) );
  SetSolverOption("output_flag", 0);
}

void HighsCommon::CloseSolver() {
  loader().Highs_destroy(lp());
}

int64_t HighsCommon::getInt64Attr(const char* name)  const {
  int64_t value = 0;
  HIGHS_CCALL(loader().Highs_getInt64InfoValue(lp(), name, &value));
  return value;
}
int HighsCommon::getIntAttr(const char* name)  const {
  int value = 0;
  HIGHS_CCALL(loader().Highs_getIntInfoValue(lp(), name, &value));
  return value;
}
double HighsCommon::getDblAttr(const char* name) const  {
  double value = 0;
  HIGHS_CCALL(loader().Highs_getDoubleInfoValue(lp(), name, &value));
  return value;
}

int HighsCommon::NumLinCons() const {
  return loader().Highs_getNumRows(lp());
}

int HighsCommon::NumVars() const {
  return loader().Highs_getNumCols(lp());
}

int HighsCommon::NumObjs()  {
  return accObjectives().numObjs();
}



void checkOption(int retvalue, const char* key) {
  if (retvalue != kHighsStatusOk)
    throw std::runtime_error(fmt::format("Error while setting option '{}'", key));
  
}
void HighsCommon::GetSolverOption(const char* key, int& value) const {
  int type;
  loader().Highs_getOptionType(lp(), key, &type);
  if (type == kHighsOptionTypeBool)
    HIGHS_CCALL(loader().Highs_getBoolOptionValue(lp(), key, &value));
  else
    HIGHS_CCALL(loader().Highs_getIntOptionValue(lp(), key, &value));
}

void HighsCommon::SetSolverOption(const char* key, int value) {
  int type;
  int ret;
  loader().Highs_getOptionType(lp(), key, &type);
  if (type == kHighsOptionTypeBool)
    ret = loader().Highs_setBoolOptionValue(lp(), key, value);
  else
    ret = loader().Highs_setIntOptionValue(lp(), key, value);
  checkOption(ret, key);
}

void HighsCommon::GetSolverOption(const char* key, double &value) const {
  HIGHS_CCALL(loader().Highs_getDoubleOptionValue(lp(), key, &value) );
}

void HighsCommon::SetSolverOption(const char* key, double value) {
  int ret = loader().Highs_setDoubleOptionValue(lp(), key, value);
  checkOption(ret, key);
}

void HighsCommon::GetSolverOption(const char* key, std::string &value) const {
  char option[256];
  HIGHS_CCALL(loader().Highs_getStringOptionValue(lp(), key, option));
  value = option;
}

void HighsCommon::SetSolverOption(const char* key, const std::string& value) {
  int ret = loader().Highs_setStringOptionValue(lp(), key, value.c_str());
  checkOption(ret, key);
}


} // namespace mp
