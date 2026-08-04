#ifndef HIGHSLOADER_H
#define HIGHSLOADER_H

#include <stdexcept>
#ifdef WIN32
#include <stdint.h>   // for standard types like uint32_t
typedef void* HMODULE;
typedef const char* LPCSTR;
typedef void* FARPROC;
extern "C" {
  HMODULE __stdcall LoadLibraryA(LPCSTR lpLibFileName);
  int __stdcall FreeLibrary(HMODULE hModule);
  FARPROC __stdcall GetProcAddress(HMODULE hModule, LPCSTR lpProcName);
}
#else
#include <dlfcn.h>
typedef void* HMODULE;
#endif

extern "C" {
// To load types and enums
#include "interfaces/highs_c_api.h"
}



namespace mp {

  // === Function pointer typedefs ===
  // Core lifecycle
  typedef void* (*Highs_create_Func)();
  typedef void   (*Highs_destroy_Func)(void*);

  // Solve
  typedef int    (*Highs_run_Func)(void*);

  // Model status / objective
  typedef int    (*Highs_getModelStatus_Func)(void*);
  typedef double (*Highs_getObjectiveValue_Func)(void*);
  typedef double (*Highs_getInfinity_Func)(void*);

  // Solutions
  typedef int    (*Highs_getSolution_Func)(void*, double*, double*, double*, double*);
  typedef int    (*Highs_setSolution_Func)(void*, const double*, const double*, const double*, const double*);

  // Basis
  typedef int    (*Highs_getBasis_Func)(void*, HighsInt*, HighsInt*);
  typedef int    (*Highs_setBasis_Func)(void*, const HighsInt*, const HighsInt*);

  // Rays
  typedef int    (*Highs_getPrimalRay_Func)(void*, HighsInt*, double*);
  typedef int    (*Highs_getDualRay_Func)(void*, HighsInt*, double*);

  // Columns (variables)
  typedef int    (*Highs_addCols_Func)(void*, HighsInt, const double*, const double*, const double*, HighsInt, const HighsInt*, const HighsInt*, const double*);
  typedef int    (*Highs_getNumCols_Func)(void*);
  typedef int    (*Highs_getColsByRange_Func)(void*, HighsInt, HighsInt, HighsInt*, double*, double*, double*, HighsInt*, HighsInt*, HighsInt*, double*);
  typedef int    (*Highs_getColsBySet_Func)(void*, HighsInt, const HighsInt*, HighsInt*, double*, double*, double*, HighsInt*, HighsInt*, HighsInt*, double*);
  typedef int    (*Highs_changeColsBoundsByMask_Func)(void*, const HighsInt*, const double*, const double*);
  typedef int    (*Highs_changeColsCostByRange_Func)(void*, HighsInt, HighsInt, const double*);
  typedef int    (*Highs_changeColsIntegralityBySet_Func)(void*, HighsInt, const HighsInt*, const HighsInt*);
  typedef int    (*Highs_passColName_Func)(void*, HighsInt, const char*);
  typedef int    (*Highs_getColIntegrality_Func)(const void*, const HighsInt, HighsInt* integrality);
  // Rows (constraints)
  typedef int    (*Highs_addRows_Func)(void*, HighsInt, const double*, const double*, HighsInt, const HighsInt*, const HighsInt*, const double*);
  typedef int    (*Highs_getNumRows_Func)(void*);

  // Quadratic / Hessian
  typedef int    (*Highs_passHessian_Func)(void*, HighsInt, HighsInt, HighsInt, const HighsInt*, const HighsInt*, const double*);

  // Multi-objectives
  typedef int    (*Highs_changeObjectiveSense_Func)(void*, int);
  typedef int    (*Highs_passLinearObjectives_Func)(void*, HighsInt, const double*, const double*, const double*, const double*, const double*, const int*);

  // Options
  typedef int    (*Highs_getOptionType_Func)(void*, const char*, int*);
  typedef int    (*Highs_getBoolOptionValue_Func)(void*, const char*, int*);
  typedef int    (*Highs_setBoolOptionValue_Func)(void*, const char*, int);
  typedef int    (*Highs_getIntOptionValue_Func)(void*, const char*, int*);
  typedef int    (*Highs_setIntOptionValue_Func)(void*, const char*, int);
  typedef int    (*Highs_getDoubleOptionValue_Func)(void*, const char*, double*);
  typedef int    (*Highs_setDoubleOptionValue_Func)(void*, const char*, double);
  typedef int    (*Highs_getStringOptionValue_Func)(void*, const char*, char*);
  typedef int    (*Highs_setStringOptionValue_Func)(void*, const char*, const char*);

  // Info / Attributes
  typedef int    (*Highs_getIntInfoValue_Func)(void*, const char*, HighsInt*);
  typedef int    (*Highs_getInt64InfoValue_Func)(void*, const char*, int64_t*);
  typedef int    (*Highs_getDoubleInfoValue_Func)(void*, const char*, double*);

  // I/O
  typedef int    (*Highs_writeModel_Func)(void*, const char*);
  typedef int    (*Highs_writeSolutionPretty_Func)(void*, const char*);
)
  typedef int    (*Highs_setCallback_Func)(void*, HighsCCallbackType, void*);
  typedef int    (*Highs_startCallback_Func)(void*, const HighsInt);
  typedef int    (*Highs_stopCallback_Func)(void*, const HighsInt);

  /**
  * Wrapper to semi-seamlessly implement dynamic library loading
  */
  class HighsLoader {
  public:
    /**
    Return the library name to be loaded (with/without CUDA, OS-dependent)
    **/
    static const char* getHighsLibraryName(bool cuda) {
      if (cuda) {
      #ifdef _WIN32
              return "highs-ampl-gpu.dll";
      #elif defined(__APPLE__)
            throw std::runtime_error("GPU-based solver not supported on MacOS");
      #else
            return "libhighs-ampl-gpu.so.1";
      #endif
      }
      else {
      #ifdef _WIN32
              return "highs-ampl.dll";
      #elif defined(__APPLE__)
              return "libhighs-ampl.1.dylib";
      #else
              return "libhighs-ampl.so.1";
      #endif
      }
    }

    HighsLoader();
    ~HighsLoader();
    /*
    * Load the specified library
    */
    bool load(const char* dllPath);
    /*
    * Unload it
    */
    void unload();

    // === HiGHS Function Pointers ===

    // Core lifecycle
    Highs_create_Func Highs_create;
    Highs_destroy_Func Highs_destroy;

    // Solve
    Highs_run_Func Highs_run;

    // Model status / objective
    Highs_getModelStatus_Func Highs_getModelStatus;
    Highs_getObjectiveValue_Func Highs_getObjectiveValue;
    Highs_getInfinity_Func Highs_getInfinity;

    // Solutions
    Highs_getSolution_Func Highs_getSolution;
    Highs_setSolution_Func Highs_setSolution;

    // Basis
    Highs_getBasis_Func Highs_getBasis;
    Highs_setBasis_Func Highs_setBasis;

    // Rays
    Highs_getPrimalRay_Func Highs_getPrimalRay;
    Highs_getDualRay_Func Highs_getDualRay;

    // Columns (variables)
    Highs_addCols_Func Highs_addCols;
    Highs_getNumCols_Func Highs_getNumCols;
    Highs_getColsByRange_Func Highs_getColsByRange;
    Highs_getColsBySet_Func Highs_getColsBySet;
    Highs_getColIntegrality_Func Highs_getColIntegrality;

    Highs_changeColsBoundsByMask_Func Highs_changeColsBoundsByMask;
    Highs_changeColsCostByRange_Func Highs_changeColsCostByRange;
    Highs_changeColsIntegralityBySet_Func Highs_changeColsIntegralityBySet;
    Highs_passColName_Func Highs_passColName;

    // Rows (constraints)
    Highs_addRows_Func Highs_addRows;
    Highs_getNumRows_Func Highs_getNumRows;

    // Quadratic / Hessian
    Highs_passHessian_Func Highs_passHessian;

    // Multi-objectives
    Highs_changeObjectiveSense_Func Highs_changeObjectiveSense;
    Highs_passLinearObjectives_Func Highs_passLinearObjectives;

    // Options
    Highs_getOptionType_Func Highs_getOptionType;
    Highs_getBoolOptionValue_Func Highs_getBoolOptionValue;
    Highs_setBoolOptionValue_Func Highs_setBoolOptionValue;
    Highs_getIntOptionValue_Func Highs_getIntOptionValue;
    Highs_setIntOptionValue_Func Highs_setIntOptionValue;
    Highs_getDoubleOptionValue_Func Highs_getDoubleOptionValue;
    Highs_setDoubleOptionValue_Func Highs_setDoubleOptionValue;
    Highs_getStringOptionValue_Func Highs_getStringOptionValue;
    Highs_setStringOptionValue_Func Highs_setStringOptionValue;

    // Info / Attributes
    Highs_getIntInfoValue_Func Highs_getIntInfoValue;
    Highs_getInt64InfoValue_Func Highs_getInt64InfoValue;
    Highs_getDoubleInfoValue_Func Highs_getDoubleInfoValue;

    // I/O
    Highs_writeModel_Func Highs_writeModel;
    Highs_writeSolutionPretty_Func Highs_writeSolutionPretty;

    // Callbacks (optional, see load())
    Highs_setCallback_Func Highs_setCallback;
    Highs_startCallback_Func Highs_startCallback;
    Highs_stopCallback_Func Highs_stopCallback;

  private:
    HMODULE hDLL;
  };

} // namespace mp





#endif // HIGHSLOADER_H