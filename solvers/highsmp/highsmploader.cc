#include "highsmploader.h"

#include "mp/format.h"

namespace mp {
HighsLoader::HighsLoader() : hDLL(nullptr) {}

HighsLoader::~HighsLoader() {
  unload();
}

bool HighsLoader::load(const char* dllPath) {
#ifdef _WIN32
  hDLL = LoadLibraryA(dllPath);
#else
  hDLL = dlopen(dllPath, RTLD_NOW);
#endif
  if (!hDLL) {
    fmt::print("Failed to load DLL: {}\n", dllPath);
    return false;
  }
  auto LoadFunc = [&](auto& funcPtr, const char* name) {
#ifdef _WIN32
    funcPtr = reinterpret_cast<std::remove_reference_t<decltype(funcPtr)>>(GetProcAddress(hDLL, name));
#else
    void* sym = dlsym(hDLL, name);
    if (!sym) {
      const char* err = dlerror();  // Grab the error message from dlsym
      fmt::print("Failed to load symbol '{}': {}\n", name, err ? err : "unknown error");
      funcPtr = nullptr;
    }
    else {
      funcPtr = reinterpret_cast<std::remove_reference_t<decltype(funcPtr)>>(sym);
    }
#endif

    if (!funcPtr)
      fmt::print("Failed to load function '{}'\n", name);
    };


  // Load the core lifecycle functions
  LoadFunc(Highs_create, "Highs_create");
  LoadFunc(Highs_destroy, "Highs_destroy");

  // Load solver functions
  LoadFunc(Highs_run, "Highs_run");
  LoadFunc(Highs_getModelStatus, "Highs_getModelStatus");
  LoadFunc(Highs_getObjectiveValue, "Highs_getObjectiveValue");
  LoadFunc(Highs_getInfinity, "Highs_getInfinity");

  // Load solutions
  LoadFunc(Highs_getSolution, "Highs_getSolution");
  LoadFunc(Highs_setSolution, "Highs_setSolution");

  // Basis
  LoadFunc(Highs_getBasis, "Highs_getBasis");
  LoadFunc(Highs_setBasis, "Highs_setBasis");

  // Rays
  LoadFunc(Highs_getPrimalRay, "Highs_getPrimalRay");
  LoadFunc(Highs_getDualRay, "Highs_getDualRay");

  // Columns
  LoadFunc(Highs_addCols, "Highs_addCols");
  LoadFunc(Highs_getNumCols, "Highs_getNumCols");
  LoadFunc(Highs_getColsByRange, "Highs_getColsByRange");
  LoadFunc(Highs_getColsBySet, "Highs_getColsBySet");
  LoadFunc(Highs_getColIntegrality, "Highs_getColIntegrality");
  LoadFunc(Highs_changeColsBoundsByMask, "Highs_changeColsBoundsByMask");
  LoadFunc(Highs_changeColsCostByRange, "Highs_changeColsCostByRange");
  LoadFunc(Highs_changeColsIntegralityBySet, "Highs_changeColsIntegralityBySet");
  LoadFunc(Highs_passColName, "Highs_passColName");

  // Rows
  LoadFunc(Highs_addRows, "Highs_addRows");
  LoadFunc(Highs_getNumRows, "Highs_getNumRows");

  // Quadratic / Hessian
  LoadFunc(Highs_passHessian, "Highs_passHessian");

  // Multi-objectives
  LoadFunc(Highs_changeObjectiveSense, "Highs_changeObjectiveSense");
  LoadFunc(Highs_passLinearObjectives, "Highs_passLinearObjectives");

  // Options
  LoadFunc(Highs_getOptionType, "Highs_getOptionType");
  LoadFunc(Highs_getBoolOptionValue, "Highs_getBoolOptionValue");
  LoadFunc(Highs_setBoolOptionValue, "Highs_setBoolOptionValue");
  LoadFunc(Highs_getIntOptionValue, "Highs_getIntOptionValue");
  LoadFunc(Highs_setIntOptionValue, "Highs_setIntOptionValue");
  LoadFunc(Highs_getDoubleOptionValue, "Highs_getDoubleOptionValue");
  LoadFunc(Highs_setDoubleOptionValue, "Highs_setDoubleOptionValue");
  LoadFunc(Highs_getStringOptionValue, "Highs_getStringOptionValue");
  LoadFunc(Highs_setStringOptionValue, "Highs_setStringOptionValue");

  // Info / Attributes
  LoadFunc(Highs_getIntInfoValue, "Highs_getIntInfoValue");
  LoadFunc(Highs_getInt64InfoValue, "Highs_getInt64InfoValue");
  LoadFunc(Highs_getDoubleInfoValue, "Highs_getDoubleInfoValue");

  // I/O
  LoadFunc(Highs_writeModel, "Highs_writeModel");
  LoadFunc(Highs_writeSolutionPretty, "Highs_writeSolutionPretty");


  auto LoadOptionalFunc = [&](auto& funcPtr, const char* name) {
#ifdef _WIN32
    funcPtr = reinterpret_cast<std::remove_reference_t<decltype(funcPtr)>>(GetProcAddress(hDLL, name));
#else
    funcPtr = reinterpret_cast<std::remove_reference_t<decltype(funcPtr)>>(dlsym(hDLL, name));
#endif
    };
  LoadOptionalFunc(Highs_setCallback, "Highs_setCallback");
  LoadOptionalFunc(Highs_startCallback, "Highs_startCallback");
  LoadOptionalFunc(Highs_stopCallback, "Highs_stopCallback");

  // Verify all pointers
  if (!(Highs_create && Highs_destroy &&
    Highs_run && Highs_getSolution && Highs_getObjectiveValue &&
    Highs_writeModel && Highs_writeSolutionPretty &&
    Highs_getIntInfoValue && Highs_getInt64InfoValue &&
    Highs_getDoubleInfoValue && Highs_getColsByRange &&
    Highs_getColsBySet && Highs_changeColsBoundsByMask &&
    Highs_setSolution && Highs_getBasis && Highs_setBasis &&
    Highs_getPrimalRay && Highs_getDualRay &&
    Highs_getColIntegrality && Highs_getModelStatus &&
    Highs_getInfinity && Highs_changeColsCostByRange &&
    Highs_changeObjectiveSense && Highs_passLinearObjectives &&
    Highs_getNumRows && Highs_getNumCols &&
    Highs_getOptionType && Highs_getBoolOptionValue &&
    Highs_getIntOptionValue && Highs_setBoolOptionValue &&
    Highs_setIntOptionValue && Highs_getDoubleOptionValue &&
    Highs_setDoubleOptionValue && Highs_getStringOptionValue &&
    Highs_setStringOptionValue && Highs_addCols &&
    Highs_changeColsIntegralityBySet && Highs_passColName &&
    Highs_passHessian && Highs_addRows)) {
    fmt::print("Failed to load one or more required functions from {}\n", dllPath);
    unload();
    return false;
  }
  return true;
}

void HighsLoader::unload() {
  if (!hDLL)
    return;
#ifdef _WIN32
  FreeLibrary(hDLL);
#else
  dlclose(hDLL);
#endif
  hDLL = nullptr;

  // Reset all function pointers to nullptr
  Highs_create = nullptr;
  Highs_destroy = nullptr;
  Highs_run = nullptr;
  Highs_getSolution = nullptr;
  Highs_getObjectiveValue = nullptr;
  Highs_writeModel = nullptr;
  Highs_writeSolutionPretty = nullptr;
  Highs_getIntInfoValue = nullptr;
  Highs_getInt64InfoValue = nullptr;
  Highs_getDoubleInfoValue = nullptr;
  Highs_getColsByRange = nullptr;
  Highs_getColsBySet = nullptr;
  Highs_changeColsBoundsByMask = nullptr;
  Highs_setSolution = nullptr;
  Highs_getBasis = nullptr;
  Highs_setBasis = nullptr;
  Highs_getPrimalRay = nullptr;
  Highs_getDualRay = nullptr;
  Highs_getColIntegrality = nullptr;
  Highs_getModelStatus = nullptr;
  Highs_getInfinity = nullptr;
  Highs_changeColsCostByRange = nullptr;
  Highs_changeObjectiveSense = nullptr;
  Highs_passLinearObjectives = nullptr;
  Highs_getNumRows = nullptr;
  Highs_getNumCols = nullptr;
  Highs_getOptionType = nullptr;
  Highs_getBoolOptionValue = nullptr;
  Highs_getIntOptionValue = nullptr;
  Highs_setBoolOptionValue = nullptr;
  Highs_setIntOptionValue = nullptr;
  Highs_getDoubleOptionValue = nullptr;
  Highs_setDoubleOptionValue = nullptr;
  Highs_getStringOptionValue = nullptr;
  Highs_setStringOptionValue = nullptr;

  Highs_addCols = nullptr;
  Highs_changeColsIntegralityBySet = nullptr;
  Highs_passColName = nullptr;
  Highs_passHessian = nullptr;
  Highs_addRows = nullptr;

  Highs_setCallback = nullptr;
  Highs_startCallback = nullptr;
  Highs_stopCallback = nullptr;
}

} // namespace mp