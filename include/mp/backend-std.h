/*
 Abstract solver backend wrapper.

 Copyright (C) 2021 AMPL Optimization Inc

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.
*/

#ifndef STD_BACKEND_H_
#define STD_BACKEND_H_

#include <chrono>
#include <cmath>
#include <functional>

#include "mp/utils-clock.h"
#include "mp/backend-with-mm.h"

/// Issue this if you redefine std feature switches
#define USING_STD_FEATURES using BaseBackend::STD_FEATURE_QUERY_FN
/// Default switch for unmentioned std features,
/// can be used to override base class' settings
#define DEFAULT_STD_FEATURES_TO( val ) \
  template <class FeatureStructName> \
  static constexpr bool STD_FEATURE_QUERY_FN( \
    const FeatureStructName& ) { return val; }
/// And default them
DEFAULT_STD_FEATURES_TO( false )
#define DEFINE_STD_FEATURE( name ) \
  struct STD_FEATURE_STRUCT_NM( name ) { };
#define ALLOW_STD_FEATURE( name, val ) \
  static constexpr bool STD_FEATURE_QUERY_FN( \
    const STD_FEATURE_STRUCT_NM( name )& ) { return val; }
#define IMPL_HAS_STD_FEATURE( name ) MP_DISPATCH_STATIC( \
  STD_FEATURE_QUERY_FN( \
    STD_FEATURE_STRUCT_NM( name )() ) )
#define STD_FEATURE_QUERY_FN AllowStdFeature__func
#define STD_FEATURE_STRUCT_NM( name ) StdFeatureDesc__ ## name

#define FEATURE_API_TO_IMPLEMENT( name, functionsignature )\
  virtual functionsignature { \
    throw std::runtime_error( \
      "Feature " #name \
        ": Function " #functionsignature " not implemented"); \
  }

// Enable driver-by-driver versioning without breaking compatibility
#ifndef DRIVER_DATE
  #define DRIVER_DATE MP_DATE
#endif



namespace mp {

/// Solution (postsolved)
struct Solution {
  /// primal
  std::vector<double> primal;
  /// dual
  std::vector<double> dual;
  /// objective values
  std::vector<double> objvals;
  /// Sparsity, if solver wants it
  ArrayRef<int> spars_primal;
};

/// StdBackend: the standard solver API wrapper
///
/// The standard wrapper provides common functionality:
/// standard option handling
/// and placeholders for solver API
template <class Impl>
class StdBackend :
    public BackendWithModelManager
{
  ////////////////////////////////////////////////////////////////////////////
  ///////////////////// TO IMPLEMENT IN THE FINAL CLASS //////////////////////
  ////////////////////////////////////////////////////////////////////////////
public:
  /// Name displayed in messages
  static const char* GetSolverName() { return "SomeSolver"; }
  /// And version
  static std::string GetSolverVersion() { return "1.0.0"; }
  /// AMPL solver name is used to parse solver options
  /// for the [name]_options environment variable.
  /// This is only done if the [executable_name]_options
  /// variable is not provided.
  static const char* GetAMPLSolverName() { return "solver"; }
  /// Unused
  static const char* GetAMPLSolverLongName() { return nullptr; }
  /// Unused
  static long Date() { return DRIVER_DATE; }


  /// Further, using virtual methods for convenience (CRTP not essential)

  /// Placeholder for solution getter (postsolved, final solution)
  virtual Solution GetSolution() = 0;
  /// Placeholder for objective values getter (postsolved, final solution)
  /// Used only if we don't want the solution values,
  /// otherwise GetSolution() provides this
  virtual ArrayRef<double> GetObjectiveValues() = 0;

  ////////////////////////////////////////////////////////////
  /////////////// OPTIONAL STANDARD FEATURES /////////////////
  /////////////// DISALLOW BY DEFAULT        /////////////////
  /////////////// IF ALLOWED, A CORR. OPTION IS CREATED //////
  /////////////// PLACEHOLDERS FOR CORR. API /////////////////
  ////////////////////////////////////////////////////////////
protected:
  /**
   * MULTIOBJ support
   */
  DEFINE_STD_FEATURE( MULTIOBJ )
  ALLOW_STD_FEATURE( MULTIOBJ, false )
  /// Placeholder: set objective priorities
  virtual void ObjPriorities(ArrayRef<int>)
  { MP_UNSUPPORTED("Backend::ObjPriorities"); }
  /// Placeholder: set objective weights.
  /// Presolve the values if needed
  virtual void ObjWeights(ArrayRef<double>) { }
  /// Placeholder: set objective abs tol
  /// Presolve the values if needed
  virtual void ObjAbsTol(ArrayRef<double>) { }
  /// Placeholder: set objective rel tol
  /// Presolve the values if needed
  virtual void ObjRelTol(ArrayRef<double>) { }


  /**
   * MULTISOL support
   * No API to overload,
   *  Impl should check need_multiple_solutions()
   *  and call ReportIntermediateSolution({x, pi, objvals}) for each
   *  (postsolve the values if needed)
   **/
  DEFINE_STD_FEATURE( MULTISOL )
  ALLOW_STD_FEATURE( MULTISOL, false )
  /**
  * Kappa estimate
  **/
  DEFINE_STD_FEATURE( KAPPA )
  ALLOW_STD_FEATURE( KAPPA, false )
  /// Placeholder: retrieve Kappa
  virtual double Kappa() { return 0.0; }
  /**
  * FeasRelax
  * No API to overload,
  * Impl should check:
  * - feasrelax() returns feasrelax mode
  * - feasrelax().(methods) give the API
  **/
  DEFINE_STD_FEATURE(FEAS_RELAX)
  ALLOW_STD_FEATURE(FEAS_RELAX, false)
  /**
  * MIP solution rounding
  * Nothing to do for the Impl, enabled by default
  **/
  DEFINE_STD_FEATURE(WANT_ROUNDING)
  ALLOW_STD_FEATURE(WANT_ROUNDING, true)


      /**
      * Model export.
      * Redefine to true and implement the method.
      **/
  DEFINE_STD_FEATURE(WRITE_PROBLEM)
  ALLOW_STD_FEATURE(WRITE_PROBLEM, false)
  FEATURE_API_TO_IMPLEMENT(WRITE_PROBLEM,
                           void DoWriteProblem(const std::string& ))
      /// Redefine this if you want some extensions
      /// to be written after solving instead.
      /// This is for compatibility wiht ASL drivers
      /// where 'writeprob' was used for native-format model
      /// and solution output #218.
      virtual std::set<std::string> NativeResultExtensions() const
  { return {".sol", ".ilp", ".mst", ".hnt", ".bas", ".json"}; }

      /**
      * Solution export.
      * Redefine to true and implement the method.
      **/
  DEFINE_STD_FEATURE(WRITE_SOLUTION)
  ALLOW_STD_FEATURE(WRITE_SOLUTION, false)
  FEATURE_API_TO_IMPLEMENT(WRITE_SOLUTION,
                           void DoWriteSolution(const std::string& ))

  ////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////// MODEL QUERY //////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  /// Helps to choose input infos
  virtual bool IsMIP() const = 0;

  ////////////////////////////////////////////////////////////////////////////
  /////////////////////////// BASIC PROCESS LOGIC ////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
public:
  /// Runs Solver given the NL file name.
  /// This is called from BackendApp.
  void RunFromNLFile(const std::string& nl_filename,
                     const std::string& filename_no_ext) override {
    ReadNL(nl_filename, filename_no_ext, GetArgvOptions());
    InputExtras();

    SetupTimerAndInterrupter();
    if (exportFileMode() > 0)
      ExportModel(export_file_names());

    // exportFileMode == 2 -> do not solve, just export
    if (exportFileMode() != 2) 
    {
      Solve();
      RecordSolveTime();

      Report();
    }
  }

  /// Detailed steps for AMPLS C API

  /// Read NL.
  /// @param opts: extra options
  /// (to be read after the [solver]_options env var).
  /// All model-related options should be here
  /// (obj:.../objno/multiobj, cvt:..., acc:...).
  void ReadNL(const std::string& nl_filename,
              const std::string& filename_no_ext,
              char** opts) override {
    GetMM().ReadNLModel(nl_filename, filename_no_ext,
                        GetCallbacks().check,
                        [this, opts](){        // options parser lambda
      this->ParseSolverOptions(opts, GetOptionFlags());
    });
  }

  /// Input warm start, suffixes, and all that can modify the model.
  /// This is used by the AMPLS C API
  void InputExtras() override {
    InputStdExtras();
  }
  void SetSolutionFileName(const std::string& filename) {
    GetMM().SetSolutionFileName(filename);
  }
  /// Report results
  void ReportResults() override {
    ReportSuffixes();
    ReportSolution();
  }


protected:
  /// Solve, no model modification any more.
  /// Can report intermediate results via ReportIntermediateSolution() during this,
  /// otherwise in ReportResults()
  virtual void Solve() override = 0;

  /// Report
  virtual void Report() {
    ReportResults();
//    if (verbose_mode())
//    if (!ampl_flag())     // otherwise to solve_message
//      PrintWarnings();
    if ( timing() )
      PrintTimingInfo();
  }

  /// Standard extras
  virtual void InputStdExtras() {
    if (multiobj()) {
      if (auto suf = ReadSuffix(suf_objpriority))
        ObjPriorities( suf );
      if (auto suf = ReadSuffix(suf_objweight))
        ObjWeights( suf );
      if (auto suf = ReadSuffix(suf_objabstol))
        ObjAbsTol( suf );
      if (auto suf = ReadSuffix(suf_objreltol))
        ObjRelTol( suf );
    }
    if (feasrelax())
      InputFeasrelaxData();
  }

  /// Timer, interrupter
  virtual void SetupTimerAndInterrupter() {
    SetupInterrupter();
    RecordSetupTime();
  }

  /// Call SetInterrupter() which can be defined in Impl
  virtual void SetupInterrupter() {
    SetInterrupter( interrupter() );
  }

  /// Placeholder for interrupt notifier in Impl
  virtual void SetInterrupter(mp::Interrupter*) = 0;

  /// Record setup time
  virtual void RecordSetupTime() {
    stats_.setup_time = GetTimeAndReset(stats_.time);
  }

  /// Record solve time
  virtual void RecordSolveTime() {
    stats_.solution_time = GetTimeAndReset(stats_.time);
  }

  /// Input feasrelax data
  /// Impl should presolve it if needed
  virtual void InputFeasrelaxData() {
    auto suf_lbpen = ReadDblSuffix( {"lbpen", suf::VAR} );
    auto suf_ubpen = ReadDblSuffix( {"ubpen", suf::VAR} );
    auto suf_rhspen = ReadDblSuffix( {"rhspen", suf::CON} );
    if (suf_lbpen.empty() && suf_ubpen.empty() && suf_rhspen.empty() &&
        0.0>lbpen() && 0.0>ubpen() && 0.0>rhspen())
      return;
    feasrelax().lbpen_ = FillFeasrelaxPenalty(suf_lbpen, lbpen(),
                GetSuffixSize(suf::VAR));
    feasrelax().ubpen_ = FillFeasrelaxPenalty(suf_ubpen, ubpen(),
                GetSuffixSize(suf::VAR));
    feasrelax().rhspen_ = FillFeasrelaxPenalty(suf_rhspen, rhspen(),
                GetSuffixSize(suf::CON));
  }

  /// For Impl to check if user wants multiple solutions reported
  using BasicSolver::need_multiple_solutions;

  /// Report suffixes
  virtual void ReportSuffixes() {
    try {
      ReportStandardSuffixes();
      ReportCustomSuffixes();
    } catch (const std::exception& exc) {
      AddWarning("SUFFIX_OUT",   // Can add warning before SOL output
                 std::string("Error reporting a suffix: ")
                     + exc.what());
    }
  }

  /// Report standard suffixes
  virtual void ReportStandardSuffixes() {
    if (IsProblemSolved() && exportKappa()) {
      ReportKappa();
    }
  }

  /// Report Kappa
  virtual void ReportKappa() {
    if (exportKappa() && 2)
    {
      double value = Kappa();
      ReportSingleSuffix(suf_objkappa, value);
      ReportSingleSuffix(suf_probkappa, value);
    }
  }

  /// Placeholder: report custom suffixes
  virtual void ReportCustomSuffixes() { }

  /// Callback.
  /// @param sol: postsolved solution.
  /// Normally report all alternative solutions,
  /// including the final one,
  /// before reporting the final solution
  /// via standard way (ReportSolution).
  void ReportIntermediateSolution(Solution sol) {
    fmt::MemoryWriter writer;
    writer.write("{}: {} {}", MP_DISPATCH( long_name() ),
                 "Alternative solution", ++kIntermSol_);
    double obj_value = std::numeric_limits<double>::quiet_NaN();
    if (sol.objvals.size()) {
      obj_value = sol.objvals[0];
      writer.write(", objective {}",
                   MP_DISPATCH( FormatObjValue(obj_value) ));
      if (obj_value > objIntermSol_.first)
        objIntermSol_.first = obj_value;
      if (obj_value < objIntermSol_.second)
        objIntermSol_.second = obj_value;
    }
    writer.write("\n");
    auto wrn = GetWarnings();
    if (wrn.size())
      writer.write("\n{}", wrn);
    const auto& scw0 = GetWarning(GetSolCheckWarningKey(false));
    const auto& scw1 = GetWarning(GetSolCheckWarningKey(true));
    if (scw0.second.size() || scw1.second.size()) {
      ++ stats_.n_altern_sol_checks_failed_;
      ClearWarning(GetSolCheckWarningKey(false));
      ClearWarning(GetSolCheckWarningKey(true));
    }
    if (round() && MPD(IsMIP()))
      RoundSolution(sol.primal, writer);
    HandleFeasibleSolution(SolveCode(), writer.c_str(),
                   sol.primal.empty() ? 0 : sol.primal.data(),
                   sol.dual.empty() ? 0 : sol.dual.data(),
                   obj_value);
  }

  /// Report final solution
  virtual void ReportSolution() {
    ReportSolution2AMPL();
    ReportSolutionViaSolver();
  }

  /// to AMPL via SOL file
  virtual void ReportSolution2AMPL() {
    double obj_value = std::numeric_limits<double>::quiet_NaN();
    auto sol = GetSolution();             // even if just dual
    fmt::MemoryWriter writer;
    writer.write("{}: {}", MP_DISPATCH( long_name() ), SolveStatus());
    if (IsProblemSolvedOrFeasible()) {
      if (sol.objvals.size()) {
        if (sol.objvals.size() > 1)
        {
          writer.write("; objective {}", FormatObjValue(sol.objvals[0]));
          writer.write("\nIndividual objective values:");
          for (size_t i = 0; i < sol.objvals.size(); i++)
            writer.write("\n\t_sobj[{}] = {}", i+1, // indexing of _sobj starts from 1
              FormatObjValue(sol.objvals[i]));
        }
        else {
          obj_value = sol.objvals[0];
          writer.write("; ");
          if (feasrelax())
            writer.write("feasrelax ");
          writer.write("objective {}",
            FormatObjValue(obj_value));
          if (feasrelax().orig_obj_available_)
            writer.write("\nOriginal objective = {}",
                         FormatObjValue(feasrelax().orig_obj_value_));
        }
      }
      if (round() && MPD(IsMIP()))
        RoundSolution(sol.primal, writer);
    }
    if (exportKappa() && 1)
      writer.write("\nkappa value: {}", Kappa());
    if (solver_msg_extra_.size()) {
      writer.write("\n");
      writer.write(solver_msg_extra_);
    }
    /// Summary on alternative solutions
    if (kIntermSol_) {
      if (objIntermSol_.first > -1e50)
        writer.write("\n{} alternative solution(s)\n"
                     "  with objective values {}..{}\n"
                     "  written to '{}1.sol' ... '{}{}.sol'.\n",
                     kIntermSol_,
                     objIntermSol_.first, objIntermSol_.second,
                     solution_stub(), solution_stub(),
                     kIntermSol_);
      else
        writer.write("\n{} alternative solution(s)\n"
                     "  written to '{}1.sol' ... '{}{}.sol'.\n",
                     kIntermSol_,
                     solution_stub(), solution_stub(),
                     kIntermSol_);
      if (stats_.n_altern_sol_checks_failed_)
        writer.write(
              "{} alternative solution checks failed.\n",
              stats_.n_altern_sol_checks_failed_);
    }
    auto wrn = GetWarnings();
    if (wrn.size())
      writer.write("\n{}", wrn);
    /// Even without a feasible solution, we can have duals/suffixes
    HandleSolution(SolveCode(), writer.c_str(),
                   sol.primal.empty() ? 0 : sol.primal.data(),
                   sol.dual.empty() ? 0 : sol.dual.data(), obj_value);
  }

  /// Abort
  virtual void Abort(int solve_code_now, std::string msg) {
    MP_RAISE_WITH_CODE(solve_code_now, msg);
  }

  /// Final timing info
  virtual void PrintTimingInfo() {
    double output_time = GetTimeAndReset(stats_.time);
    MP_DISPATCH( Print("Setup time = {:.6f}s\n"
                       "Solution time = {:.6f}s\n"
                       "Output time = {:.6f}s\n",
                       stats_.setup_time, stats_.solution_time, output_time) );
  }

  /// MIP solution rounding. Better use solver's capabilities
  void RoundSolution(std::vector<double>& sol,
                     fmt::MemoryWriter& writer) {
    auto rndres = DoRound(sol);
    if (rndres.first) {
      ModifySolveCodeAndMessageAfterRounding(rndres, writer);
    }
  }

  std::pair<int, double> DoRound(std::vector<double>& sol) {
    int nround = 0;
    double maxmodif = 0.0;
    const bool fAssign = round() & 1;
    const auto& fInt = IsVarInt();
    for (auto j = std::min(fInt.size(), sol.size()); j--; ) {
      if (fInt[j]) {
        auto y = std::round(sol[j]);
        auto d = std::fabs(sol[j]-y);
        if (0.0!=d) {
          ++nround;
          if (maxmodif < d)
            maxmodif = d;
          if (fAssign)
            sol[j] = y;
        }
      }
    }
    return {nround, maxmodif};
  }

  void ModifySolveCodeAndMessageAfterRounding(
        std::pair<int, double> rndres, fmt::MemoryWriter& writer) {
    if (round() & 2 && IsSolStatusRetrieved()) {
      // old: solve_code_ = 3 - (round() & 1);
    }
    if (round() & 4) {
      auto sc = rndres.first > 1 ? "s" : "";
      writer.write(
            "\n{} integer variable{} {}rounded to integer{};"
      " maxerr = {:.16}", rndres.first, sc,
            round() & 1 ? "" : "would be ", sc, rndres.second);
    }
  }

  //////////////////////// SOLUTION STATUS ACCESS ///////////////////////////////

  /// Solve result number
  virtual int SolveCode() const { return status_.first; }
  /// Solver result message
  const char* SolveStatus() const { return status_.second.c_str(); }
  /// Set solve result
  void SetStatus(std::pair<int, std::string> stt) { status_=stt; }

  //////////////////////// SOLUTION STATUS ADAPTERS ///////////////////////////////
  /** Following the taxonomy of the enum sol, returns true if
      we have an optimal solution or a feasible solution for a
      satisfaction problem */
  virtual bool IsProblemSolved() const {
    assert(IsSolStatusRetrieved());
    return sol::SOLVED<=SolveCode()
        && SolveCode()<=sol::SOLVED_LAST;
  }
  /// Solved or feasible
  virtual bool IsProblemSolvedOrFeasible() const {
    assert( IsSolStatusRetrieved() );
    return
        (sol::SOLVED_LAST>=SolveCode() &&
         sol::SOLVED<=SolveCode())
        ||
        (sol::LIMIT_FEAS>=SolveCode() &&
         sol::LIMIT_FEAS_LAST<=SolveCode())
        ||
        (sol::UNBOUNDED_FEAS>=SolveCode() &&
         sol::UNBOUNDED_NO_FEAS_LAST<=SolveCode());
  }
  /// Undecidedly infeas or unbnd
  virtual bool IsProblemIndiffInfOrUnb() const {
    assert( IsSolStatusRetrieved() );
    return sol::LIMIT_INF_UNB<=SolveCode()
        && SolveCode()<=sol::LIMIT_INF_UNB_LAST;
  }
  /// Infeasible or unbounded
  virtual bool IsProblemInfOrUnb() const {
    assert( IsSolStatusRetrieved() );
    auto sc = SolveCode();
    return
        (sol::INFEASIBLE<=sc
         && sol::UNBOUNDED_NO_FEAS_LAST>=sc)
        || IsProblemIndiffInfOrUnb();
  }
  virtual bool IsProblemInfeasible() const {
    assert( IsSolStatusRetrieved() );
    auto sc = SolveCode();
    return sol::INFEASIBLE<=sc && sol::INFEASIBLE_LAST>sc;
  }
  virtual bool IsProblemUnbounded() const {
    assert( IsSolStatusRetrieved() );
    auto sc = SolveCode();
    return sol::UNBOUNDED_FEAS<=sc
        && sol::UNBOUNDED_NO_FEAS_LAST>=sc;
  }
  virtual bool IsSolStatusRetrieved() const {
    return sol::NOT_SET!=SolveCode();
  }

  struct Stats {
    std::chrono::steady_clock::time_point time = std::chrono::steady_clock::now();
    double setup_time = 0.0;
    double solution_time = 0.0;
    int n_altern_sol_checks_failed_ = 0;
  };
  Stats stats_;


  /////////////////////////////// SOME MATHS ////////////////////////////////

  /// Use the solver' inf
  bool IsFinite(double n) const {
    return n>MP_DISPATCH( MinusInfinity() ) &&
        n<MP_DISPATCH( Infinity() );
  }

  /// AMPL's inf
  static constexpr double AMPLInf()
  { return INFINITY; }
  /// AMPL's -inf
  static constexpr double AMPLMinusInf()
  { return -INFINITY; }


public:
  using BasicSolver::Print;
  using BasicSolver::add_to_long_name;
  using BasicSolver::add_to_version;
  using BasicSolver::set_option_header;
  using BasicSolver::add_to_option_header;


  ///////////////////////// STORING SOLUTON STATUS //////////////////////
private:
  std::pair<int, std::string> status_ { sol::NOT_SET, "status not set" };
  int kIntermSol_ = 0;   // last written intermed solution index
  std::pair<double, double> objIntermSol_ { -1e100, 1e100 };

  ///////////////////////// STORING SOLVER MESSAGES //////////////////////
private:
  std::string solver_msg_extra_;

protected:
  void AddToSolverMessage(const std::string& msg)
  { solver_msg_extra_ += msg; }

  ///////////////////////////// OPTIONS /////////////////////////////////
public:
  using BasicSolver::AddOption;
  using BasicSolver::AddOptionSynonyms_Inline_Front;
  using BasicSolver::AddOptionSynonyms_Inline_Back;
  using BasicSolver::AddOptionSynonyms_OutOfLine;
  using BasicSolver::FindOption;

private:
  /// Recorded solver options
  using SlvOptionRecord = std::function<void(void)>;
  std::vector< SlvOptionRecord > slvOptionRecords_;

protected:
  // set to true in constructor to disable setting options 
  // while parsing them. They can be set later with ReplaySolverOptions
  bool onlyRecordOptions_ = false; 
  // Accessed from SolverOptionAccessor to decide if to skip
  // setting options
  bool onlyRecordOptions() { return onlyRecordOptions_; }

  void RecordSolverOption(SlvOptionRecord sor)
  { slvOptionRecords_.push_back(sor); }
  void ReplaySolverOptions() {
    for (auto f: slvOptionRecords_)
      f();
  }

  /// Solver options accessor, facilitates calling
  /// backend_.Get/SetSolverOption()
  template <class Value, class Index>
  class SolverOptionAccessor {
    Impl& backend_;
  public:
    using value_type = Value;
    using index_type = Index;
    SolverOptionAccessor(Impl& b) : backend_(b) { }
    /// Options setup
    Value get(const SolverOption& , Index i) const {
      Value v;
      backend_.GetSolverOption(i, v);
      return v;
    }
    void set(const SolverOption& ,
             typename internal::OptionHelper<Value>::Arg v,
             Index i) {
      auto *pBackend = &backend_;
      auto setter = [=]() { pBackend->SetSolverOption(i, v); };
      setter();    // run it now
      backend_.RecordSolverOption(setter);
    }
  };

  template <class ValueType, class KeyType>
  class ConcreteOptionWrapper :
      public SolverOptionManager::ConcreteOptionWithInfo<
      SolverOptionAccessor<ValueType, KeyType>, ValueType, KeyType> {

    using COType = SolverOptionManager::ConcreteOptionWithInfo<
    SolverOptionAccessor<ValueType, KeyType>, ValueType, KeyType>;
    using SOAType = SolverOptionAccessor<ValueType, KeyType>;

    SOAType soa_;
  public:
    ConcreteOptionWrapper(Impl* impl_, const char *name, const char *description,
                          KeyType k, ValueArrayRef values = ValueArrayRef()) :
      COType(name, description, &soa_, &SOAType::get, &SOAType::set, k, values),
      soa_(*impl_)
    { }
  };

  /// Adding solver options of types int/double/string/...
  /// The type is deduced from the two last parameters min, max
  /// (currently unused otherwise.)
  /// If min/max omitted, assume ValueType=std::string
  /// Assumes existence of Impl::Get/SetSolverOption(KeyType, ValueType(&))
  template <class KeyType, class ValueType=std::string>
  void AddSolverOption(const char *name_list, const char *description,
                       KeyType k,
                       ValueType , ValueType ) {
    AddOption(SolverOptionManager::OptionPtr(
                new ConcreteOptionWrapper<
                ValueType, KeyType>(
                  (Impl*)this, name_list, description, k)));
  }

  /// String-valued option
  template <class KeyType, class ValueType = std::string>
  void AddSolverOption(const char* name_list, const char* description,
    KeyType k) {
    AddOption(SolverOptionManager::OptionPtr(
      new ConcreteOptionWrapper<
      ValueType, KeyType>(
        (Impl*)this, name_list, description, k)));
  }

  /// With value table.
  /// Type deduced from \a defaultValue
  template <class KeyType, class ValueType = std::string>
  void AddSolverOption(const char* name_list, const char* description,
      KeyType k, ValueArrayRef values, ValueType defaultValue) {
    internal::Unused(defaultValue);
    AddOption(SolverOptionManager::OptionPtr(
      new ConcreteOptionWrapper<
      ValueType, KeyType>(
        (Impl*)this, name_list, description, k, values)));
  }


  /// With value table: specialize for const char*.
  template <class KeyType>
  void AddSolverOption(const char* name_list, const char* description,
                       KeyType k, ValueArrayRef values,
                       const char* defaultValue) {
    internal::Unused(defaultValue);
    AddOption(SolverOptionManager::OptionPtr(
                new ConcreteOptionWrapper<
                std::string, KeyType>(
                  (Impl*)this, name_list, description, k, values)));
  }


private:
  struct Options {
    int exportKappa_ = 0;

    /// feasrelax penalty options
    double lbpen_=1.0, ubpen_=1.0, rhspen_=1.0;

    int round_=0;
    double round_reptol_=1e-9;

    /// For write prob
    std::vector<std::string> export_files_;
    std::vector<std::string> just_export_files_;
    /// For write sol
    std::vector<std::string> export_sol_files_;
  } storedOptions_;


  /// Once Impl allows FEASRELAX,
  /// it should check this via feasrelax()
  class FeasrelaxIO {
  public:
    /// Mode!=1: if & how feasrelax should be done
    operator int() const { return mode_; }
    /// Penalty vectors
    ArrayRef<double> lbpen() const { return lbpen_; }
    ArrayRef<double> ubpen() const { return ubpen_; }
    ArrayRef<double> rhspen() const { return rhspen_; }
    /// Call me if orig_obj_value() will be set
    void flag_orig_obj_available() { orig_obj_available_=true; }
    /// Access original objective value
    double& orig_obj_value() { return orig_obj_value_; }
    friend class StdBackend;
  private:
    /// --------------------- INPUT ----------------------
    int mode_=0;  // whether Impl should do it and which mode,
                  // can be redefined by Impl if cannot map
                  // from standard options
    /// Empty vector means +inf penalties
    std::vector<double> lbpen_, ubpen_, rhspen_;

    /// --------------------- OUTPUT, filled by Impl -----
    bool orig_obj_available_ = false;
    double orig_obj_value_ = 0.0;
  } feasRelaxIO_;


protected:  //////////// Option accessors ////////////////
  int exportKappa() const { return storedOptions_.exportKappa_; }

  /// Feasrelax I/O data
  FeasrelaxIO& feasrelax() { return feasRelaxIO_; }
  double lbpen() const { return storedOptions_.lbpen_; }
  double ubpen() const { return storedOptions_.ubpen_; }
  double rhspen() const { return storedOptions_.rhspen_; }

  /// Whether to round MIP solution and modify messages
  int round() const
  { return IMPL_HAS_STD_FEATURE(WANT_ROUNDING) ? storedOptions_.round_ : 0; }
  /// MIP solution rounding reporting tolerance
  double round_reptol() const { return storedOptions_.round_reptol_; }

  /// Return 2 if we only need to write the problem, 1 if solution 
  /// and export are wanted, 0 if no export is needed or supported
  int exportFileMode() const
  {
    if (IMPL_HAS_STD_FEATURE(WRITE_PROBLEM)) {
      if (!storedOptions_.export_files_.empty())
        return 1;
      if (!storedOptions_.just_export_files_.empty())
        return 2;
    }
    return 0;
  }
  const std::vector<std::string>& export_file_names() const {
    return
        storedOptions_.export_files_.empty()
        ? storedOptions_.just_export_files_
        : storedOptions_.export_files_;
  }


protected:
  virtual void InitStandardOptions() {
    if (IMPL_HAS_STD_FEATURE(KAPPA))
      AddStoredOption("alg:kappa kappa basis_cond",
        "Whether to return the estimated condition number (kappa) of "
        "the optimal basis (default 0): sum of 1 = report kappa in the result message; "
        "2 = return kappa in the solver-defined suffix .kappa on the objective and "
        "problem. The request is ignored when there is no optimal basis.",
        storedOptions_.exportKappa_);

    if (IMPL_HAS_STD_FEATURE(FEAS_RELAX)) {
      AddStoredOption("alg:feasrelax feasrelax",
                      "Whether to modify the problem into a feasibility "
                          "relaxation problem:\n"
                          "\n"
                          "| 0 = No (default)\n"
                          "| 1 = Yes, minimizing the weighted sum of violations\n"
                          "| 2 = Yes, minimizing the weighted sum of squared violations\n"
                          "| 3 = Yes, minimizing the weighted count of violations\n"
                          "| 4-6 = Same objective as 1-3, but also optimize the "
                             "original objective, subject to the violation "
                             "objective being minimized.\n"
                      "\n"
                      "Weights are given by suffixes .lbpen and .ubpen on variables "
                      "and .rhspen on constraints (when nonnegative, default values = 0), "
                      "else by keywords alg:lbpen, alg:ubpen, and alg:rhspen, respectively "
                      "(default values = 1). Weights < 0 are treated as Infinity, allowing "
                      "no violation.",
          feasrelax().mode_);
      AddStoredOption("alg:lbpen lbpen", "See alg:feasrelax.",
          storedOptions_.lbpen_);
      AddStoredOption("alg:ubpen ubpen", "See alg:feasrelax.",
          storedOptions_.ubpen_);
      AddStoredOption("alg:rhspen rhspen", "See alg:feasrelax.",
          storedOptions_.rhspen_);
    }

    if (IMPL_HAS_STD_FEATURE( WANT_ROUNDING )) {
      AddStoredOption("mip:round round",
                      "Whether to round integer variables to integral values before "
                      "returning the solution, and whether to report that the solver "
                      "returned noninteger values for integer values:  sum of\n"
                      "\n"
                      "|  1 ==> Round nonintegral integer variables\n"
                      "|  2 ==> Modify solve_result\n"
                      "|  4 ==> Modify solve_message\n"
                      "\n"
                      "Default = 0.  Modifications that were or would be made are "
                      "reported in solve_result and solve_message only if the maximum "
                      "deviation from integrality exceeded mip:round_reptol.",
                    storedOptions_.round_);
      AddStoredOption("mip:round_reptol round_reptol",
                      "Tolerance for reporting rounding of integer variables to "
                      "integer values; see \"mip:round\".  Default = 1e-9.",
                    storedOptions_.round_reptol_);
    }

    if (IMPL_HAS_STD_FEATURE(WRITE_PROBLEM)) {
      AddListOption("tech:writemodel writeprob writemodel tech:exportfile",
        "Specifies files where to export the model before "
        "solving (repeat the option for several files.) "
        "File name extensions can be ``.lp[.7z]``, ``.mps``, etc.",
        storedOptions_.export_files_);

      AddListOption("tech:writemodelonly justwriteprob justwritemodel",
        "Specifies files where to export the model, no solving "
                    "(option can be repeated.) "
        "File extensions can be ``.dlp``, ``.mps``, etc.",
        storedOptions_.just_export_files_);
    }

    if (IMPL_HAS_STD_FEATURE(WRITE_SOLUTION))
      AddListOption("tech:writesolution writesol writesolution",
        "Specifies the names of files where to export the solution "
        "and/or other result files in solver's native formats. "
        "Option can be repeated. "
        "File name extensions can be "
        "``.sol[.tar.gz]``, ``.json``, ``.bas``, ``.ilp``, etc.",
        storedOptions_.export_sol_files_);
  }

  virtual void InitCustomOptions() { }


  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////// STANDARD SUFFIXES ///////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
private:
  const SuffixDef<int> suf_objpriority = { "objpriority", suf::OBJ | suf::INPUT };
  const SuffixDef<double> suf_objweight = { "objweight", suf::OBJ | suf::INPUT };
  const SuffixDef<double> suf_objabstol = { "objabstol", suf::OBJ | suf::INPUT };
  const SuffixDef<double> suf_objreltol = { "objreltol", suf::OBJ | suf::INPUT };

  const SuffixDef<double> suf_objkappa = { "kappa", suf::OBJ | suf::OUTONLY };
  const SuffixDef<double> suf_probkappa = { "kappa", suf::PROBLEM | suf::OUTONLY };


  /////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////// SERVICE STUFF ///////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////
public:
  /// Default mp::Solver flags,
  /// used there to implement multiobj and .nsol
  static int Flags() {
    int flg=0;
    if ( IMPL_HAS_STD_FEATURE(MULTISOL) )
      flg |= BasicSolver::MULTIPLE_SOL;
    if ( IMPL_HAS_STD_FEATURE(MULTIOBJ) )
      flg |= BasicSolver::MULTIPLE_OBJ;
    return flg;
  }

  void InitMetaInfoAndOptions() override {
    BasicSolver::InitMetaInfoAndOptions(MP_DISPATCH( GetAMPLSolverName() ),
                                        "",             // Filled below
                           Date(), Flags());
    InitNamesAndVersion();
    InitStandardOptions();
    InitCustomOptions();
  }

  virtual void InitNamesAndVersion() {
    auto name = MP_DISPATCH( GetSolverName() );
    auto version = MP_DISPATCH( GetSolverVersion() );
    this->set_long_name( fmt::format("{} {}", name, version ) );
    this->set_version( fmt::format("AMPL/{} Optimizer [{}]",
                                   name, version ) );
    auto lic_cb = GetCallbacks().additionalText;
    if (lic_cb)
      this->set_license_info(lic_cb());
  }

  using MPSolverBase = BackendWithModelManager;


public:
  using MPUtils = MPSolverBase;              // Allow Converter access the SolverImpl
  const MPUtils& GetMPUtils() const { return *this; }
  MPUtils& GetMPUtils() { return *this; }

  using MPSolverBase::debug_mode;


protected:
  using MPSolverBase::interrupter;


protected:
  /// Returns {} if these penalties are +inf
  std::vector<double> FillFeasrelaxPenalty(
      ArrayRef<double> suf_pen, double pen, int n) {
    if (suf_pen.empty() && pen<0.0)
      return {};
    std::vector<double> result(n,
                               pen<0.0 ? MP_DISPATCH(Infinity()) : pen);
    for (size_t i=suf_pen.size(); i--; ) {
      result[i] = suf_pen[i]<0.0 ? MP_DISPATCH(Infinity()) : suf_pen[i];
    }
    return result;
  }


public:
  /// Write solution result in solver native format
  /// via solver's output.
  virtual void ReportSolutionViaSolver() {
    if (IMPL_HAS_STD_FEATURE(WRITE_SOLUTION))
      if (!storedOptions_.export_sol_files_.empty())
        for (const auto& fln: storedOptions_.export_sol_files_)
          try {
            DoWriteSolution(fln);
          } catch (const std::exception& exc) {
            // Not AddWarning() any more
            // because this is after solution report
            Print("    --- Error saving '{}': {}\n", fln, exc.what());
          }
  }

  /// Write model
  virtual void ExportModel(
        const std::vector<std::string>& filenames) {
    auto resext = NativeResultExtensions();
    for (const auto& fln: filenames) {
      bool fPostpone = false;
      if (IMPL_HAS_STD_FEATURE(WRITE_SOLUTION)) {
        auto dot = fln.rfind('.');  // Postpone this file?
        if (std::string::npos != dot) {
          auto ext = fln.substr(dot, fln.size()-dot);
          if (resext.end() != resext.find(ext)) {
            fPostpone = true;
          }
        }
      }
      if (fPostpone) {   // ASL-style usage of 'writeprob' #218
        storedOptions_.export_sol_files_.push_back(fln);
      } else {
        try {
          DoWriteProblem(fln);
        } catch (const std::exception& exc) {
          auto msg
              = std::string("Model export to file '")
              + fln
              + "' failed:\n"
              + exc.what();
          if (IMPL_HAS_STD_FEATURE(WRITE_SOLUTION))
            msg +=
                "\n    Note: to export solutions and results\n"
                "    in the solver's native formats,\n"
                "    use option 'tech:writesolution'";
          MP_RAISE(msg);
        }
      }
    }
  }

  /// Virtual destructor
  virtual ~StdBackend() { }
};

}  // namespace mp

#endif  // STD_BACKEND_H_
