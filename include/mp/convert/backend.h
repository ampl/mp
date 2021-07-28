/*
 Abstract solver backend wrapper.

 Copyright (C) 2020 AMPL Optimization Inc

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

 Author: Gleb Belov <gleb.belov@monash.edu>
 */

#ifndef BACKEND_H_
#define BACKEND_H_

#include <stdexcept>

#include "mp/clock.h"
#include "mp/convert/converter_query.h"
#include "mp/convert/constraint_keeper.h"
#include "mp/convert/std_constr.h"
#include "mp/convert/std_obj.h"
#include "mp/convert/model.h"
#include "mp/convert/model_adapter.h"

#define DEFINE_STD_FEATURE( name, defval ) \
  struct STD_FEATURE_STRUCT_NM( name ) { }; \
  ALLOW_STD_FEATURE( name, defval )
#define ALLOW_STD_FEATURE( name, val ) \
  static constexpr bool STD_FEATURE_QUERY_FN( \
    const STD_FEATURE_STRUCT_NM( name )& ) { return val; }
#define IMPL_HAS_STD_FEATURE( name ) MP_DISPATCH( \
  STD_FEATURE_QUERY_FN( STD_FEATURE_STRUCT_NM( name )() ) )
#define STD_FEATURE_QUERY_FN AllowStdFeature__func
#define STD_FEATURE_STRUCT_NM( name ) StdFeatureDesc__ ## name

#define RAISE_NOT_IMPLEMENTED(name) \
  throw std::runtime_error( #name  " has not been implemented!")

namespace mp {

/// Basic backend wrapper.
/// The basic wrapper provides common functionality: option handling
/// and placeholders for solver API
template <class Impl>
class BasicBackend :
    public BasicConstraintAdder,
    private SolverImpl< ModelAdapter< BasicModel<> > >   // mp::Solver stuff, hidden
{
  ConverterQuery *p_converter_query_object = nullptr;
  using MPSolverBase = SolverImpl< ModelAdapter< BasicModel<> > >;
public:
  using MPUtils = MPSolverBase;              // Allow Converter access the SolverImpl
  const MPUtils& GetMPUtils() const { return *this; }
  MPUtils& GetMPUtils() { return *this; }
public:
  BasicBackend() :
    MPSolverBase(
      Impl::GetSolverInvocationName(),
      Impl::GetAMPLSolverLongName(),
      Impl::Date(), Impl::Flags())
  { }
  virtual ~BasicBackend() { }

  void OpenSolver() { }
  void CloseSolver() { }

  /// Converter should provide this before Backend can run solving
  void ProvideConverterQueryObject(ConverterQuery* pCQ) { p_converter_query_object = pCQ; }

private: // hiding this detail, it's not for the final backends
  const ConverterQuery& GetCQ() const {
    assert(nullptr!=p_converter_query_object);
    return *p_converter_query_object;
  }
  ConverterQuery& GetCQ() {
    assert(nullptr!=p_converter_query_object);
    return *p_converter_query_object;
  }

public:

  /// Default metadata
  static const char* GetSolverName() { return "SomeSolver"; }
  static std::string GetSolverVersion() { return "-285.68.53"; }
  static const char* GetSolverInvocationName() { return "solverdirect"; }
  static const char* GetAMPLSolverLongName() { return nullptr; }
  static const char* GetBackendName()    { return "BasicBackend"; }
  static const char* GetBackendLongName() { return nullptr; }
  static long Date() { return MP_DATE; }

  /// Default flags
  static int Flags() {
    int flg=0;
    if (Impl::IfMultipleSol() )
      flg |= Solver::MULTIPLE_SOL;
    if (Impl::IfMultipleObj() )
      flg |= Solver::MULTIPLE_OBJ;
    return flg;
  }
  static bool IfMultipleSol() { return false; }
  static bool IfMultipleObj() { return false; }

  void InitMetaInfoAndOptions() {
    MP_DISPATCH( InitNamesAndVersion() );
    MP_DISPATCH( InitStandardOptions() );
    MP_DISPATCH( InitCustomOptions() );
  }

  void InitNamesAndVersion() {
    auto name = MP_DISPATCH( GetSolverName() );
    auto version = MP_DISPATCH( GetSolverVersion() );
    this->set_long_name( fmt::format("{} {}", name, version ) );
    this->set_version( fmt::format("AMPL/{} Optimizer [{}]",
                                   name, version ) );
  }

  ///////////////////////////// MODEL MANIP //////////////////////////////
  using Model = BasicModel<>;

  using Variable = typename Model::Variable;

  void InitProblemModificationPhase() { }
  void FinishProblemModificationPhase() { }
  void AddVariable(Variable var) {
    throw MakeUnsupportedError("BasicBackend::AddVariable");
  }
  void AddCommonExpression(Problem::CommonExpr cexpr) {
    throw MakeUnsupportedError("BasicBackend::AddCommonExpressions");
  }
  void AddLogicalConstraint(Problem::LogicalCon lcon) {
    throw MakeUnsupportedError("BasicBackend::AddLogicalConstraints");
  }

  void AddObjective(typename Model::Objective obj) {
    if (obj.nonlinear_expr()) {
      MP_DISPATCH( AddGeneralObjective( obj ) );
    } else {
      LinearExprUnzipper leu(obj.linear_expr());
      LinearObjective lo { obj.type(),
            std::move(leu.c_), std::move(leu.v_) };
      if (nullptr==obj.p_extra_info()) {
        MP_DISPATCH( SetLinearObjective( obj.index(), lo ) );
      } else {
        auto qt = obj.p_extra_info()->qt_;
        assert(!qt.empty());
        MP_DISPATCH( SetQuadraticObjective( obj.index(),
                       QuadraticObjective{std::move(lo), std::move(qt)} ) );
      }
    }
  }
  void AddGeneralObjective(typename Model::Objective ) {
    throw MakeUnsupportedError("BasicBackend::AddGeneralObjective");
  }
  void SetLinearObjective( int, const LinearObjective& ) {
    throw MakeUnsupportedError("BasicBackend::AddLinearObjective");
  }
  void SetQuadraticObjective( int, const QuadraticObjective& ) {
    throw MakeUnsupportedError("BasicBackend::AddQuadraticObjective");
  }

  void AddAlgebraicConstraint(typename Model::AlgebraicCon con) {
    if (con.nonlinear_expr()) {
      MP_DISPATCH( AddGeneralConstraint( con ) );
    } else {
      LinearExprUnzipper leu(con.linear_expr());
      auto lc = LinearConstraint{
          std::move(leu.c_), std::move(leu.v_),
          con.lb(), con.ub() };
      if (nullptr==con.p_extra_info()) {
        MP_DISPATCH( AddConstraint( lc ) );
      } else {
        auto qt = con.p_extra_info()->qt_;
        assert(!qt.empty());
        MP_DISPATCH( AddConstraint( QuadraticConstraint{std::move(lc), std::move(qt)} ) );
      }
    }
  }

  void AddGeneralConstraint(typename Model::AlgebraicCon ) {
    throw MakeUnsupportedError("BasicBackend::AddGeneralConstraint");
  }

  ////////////////// Some basic custom constraints /////////////////
  USE_BASE_CONSTRAINT_HANDLERS(BasicConstraintAdder)

  /// Optionally exclude LDCs from being posted,
  /// then all those are converted to LinearConstraint's first
  ACCEPT_CONSTRAINT(LinearDefiningConstraint, NotAccepted)
  void AddConstraint(const LinearDefiningConstraint& ldc) {
    MP_DISPATCH( AddConstraint(ldc.to_linear_constraint()) );
  }

  ACCEPT_CONSTRAINT(LinearConstraint, Recommended)
  /// TODO Attributes (lazy/user cut, etc)
  void AddConstraint(const LinearConstraint& ) {
    throw MakeUnsupportedError("BasicBackend::AddLinearConstraint");
  }


  ////////////////////////////////////////////////////////////////////////////
  /////////////////////////// BASIC PROCESS LOGIC ////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  void SolveAndReport() {
    MP_DISPATCH( ReadSuffixes() );

    MP_DISPATCH( PrepareSolve() );
    MP_DISPATCH( SolveAndReportIntermediateResults() );
    MP_DISPATCH( WrapupSolve() );

    MP_DISPATCH( ObtainSolutionStatus() );
    MP_DISPATCH( CalculateAndReportDerivedResults() );
    MP_DISPATCH( ReportSolution() );
    if (MP_DISPATCH( timing() ))
      MP_DISPATCH( PrintTimingInfo() );
  }

  void ReadSuffixes() {
    MP_DISPATCH( ReadStandardSuffixes() );
    MP_DISPATCH( ReadCustomSuffixes() );
  }

  void ReadStandardSuffixes() {
    if (storedOptions_.importPriorities_)
      MP_DISPATCH( VarPriorities( ReadSuffix(suf_varpriority) ) );
    if (multiobj()) {
      MP_DISPATCH( ObjPriorities( ReadSuffix(suf_objpriority) ) );
      MP_DISPATCH( ObjWeights( ReadSuffix(suf_objweight) ) );
      MP_DISPATCH( ObjAbsTol( ReadSuffix(suf_objabstol) ) );
      MP_DISPATCH( ObjRelTol( ReadSuffix(suf_objreltol) ) );
    }
  }
  void ReadCustomSuffixes() { }

  void PrepareSolve() {
    MP_DISPATCH( SetupInterrupter() );
    MP_DISPATCH( SetupTimer() );
  }

  void SetupInterrupter() {
    MP_DISPATCH( SetInterrupter(MP_DISPATCH( interrupter() )) );
  }

  void SetupTimer() {
    stats.setup_time = GetTimeAndReset(stats.time);
  }

  void WrapupSolve() {
    stats.solution_time = GetTimeAndReset(stats.time);
  }

  void ObtainSolutionStatus() {
    solve_status = MP_DISPATCH(
          ConvertSolutionStatus(*MP_DISPATCH( interrupter() ), solve_code) );
  }
  void CalculateAndReportDerivedResults() { }

  using Solver::need_multiple_solutions;

  void ReportSolution() {
    MP_DISPATCH( ReportSuffixes() );
    if (need_multiple_solutions())
      MP_DISPATCH( ReportMultipleSolutions() );
    MP_DISPATCH( ReportPrimalDualValues() );
  }

  void ReportSuffixes() {
    MP_DISPATCH( ReportStandardSuffixes() );
    MP_DISPATCH( ReportCustomSuffixes() );
  }

  void ReportStandardSuffixes() { }

  void ReportCustomSuffixes() { }

  void ReportMultipleSolutions() {
    MP_DISPATCH( StartPoolSolutions() );
    while (MP_DISPATCH( SelectNextPoolSolution() )) {
      MP_DISPATCH( ReportCurrentPoolSolution() );
    }
    MP_DISPATCH( EndPoolSolutions() );
  }

  /// When IfMultipleSol() returns true, Impl has to define these
  void StartPoolSolutions() { RAISE_NOT_IMPLEMENTED("StartPoolSolutions()"); }
  bool SelectNextPoolSolution() { return false; }   // none by default
  void EndPoolSolutions() { RAISE_NOT_IMPLEMENTED("EndPoolSolutions()"); }
  std::vector<double> CurrentPoolPrimalSolution()
  { RAISE_NOT_IMPLEMENTED("CurrentPoolPrimalSolution()"); return {}; }
  double CurrentPoolObjectiveValue() const
  { RAISE_NOT_IMPLEMENTED("CurrentPoolObjectiveValue()"); return 0.0; }

  void ReportCurrentPoolSolution() {
    double obj_value = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> solution;

    fmt::MemoryWriter writer;
    writer.write("{}: {}", MP_DISPATCH( long_name() ), solve_status);
    if (MP_DISPATCH( NumberOfObjectives() ) > 0) {
      obj_value = MP_DISPATCH( CurrentPoolObjectiveValue() );
      writer.write("; objective {}",
                   MP_DISPATCH( FormatObjValue(obj_value) ));
    }
    writer.write("\n");

    solution = MP_DISPATCH( CurrentPoolPrimalSolution() );
    HandleFeasibleSolution(writer.c_str(),
                   solution.empty() ? 0 : solution.data(),
                   nullptr, obj_value);
  }

  void ReportPrimalDualValues() {
    double obj_value = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> solution, dual_solution;

    fmt::MemoryWriter writer;
    writer.write("{}: {}", MP_DISPATCH( long_name() ), solve_status);
    if (solve_code < sol::INFEASIBLE) {
      if (MP_DISPATCH( NumberOfObjectives() ) > 0) {
        obj_value = MP_DISPATCH( ObjectiveValue() );
        writer.write("; objective {}",
                     MP_DISPATCH( FormatObjValue(obj_value) ));
      }
    }
    writer.write("\n");

    solution = MP_DISPATCH( PrimalSolution() );
    dual_solution = MP_DISPATCH( DualSolution() );  // Try in any case
    HandleSolution(solve_code, writer.c_str(),
                   solution.empty() ? 0 : solution.data(),
                   dual_solution.empty() ? 0 : dual_solution.data(), obj_value);
  }

  /// Dual solution. Returns empty if not available
  std::vector<double> DualSolution() { return {}; }

  void PrintTimingInfo() {
    double output_time = GetTimeAndReset(stats.time);
    MP_DISPATCH( Print("Setup time = {:.6f}s\n"
                       "Solution time = {:.6f}s\n"
                       "Output time = {:.6f}s\n",
                       stats.setup_time, stats.solution_time, output_time) );
  }

  /////////////////////////////// SERVICE STUFF ///////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////// SOLUTION STATUS /////////////////////////////////
  bool IsProblemInfOrUnb() const {
    assert( IsSolStatusRetrieved() );
    return sol::INFEASIBLE<=solve_code &&
        sol::UNBOUNDED>=solve_code;
  }

  bool IsProblemInfeasible() const {
    assert( IsSolStatusRetrieved() );
    return sol::INFEASIBLE<=solve_code &&
        sol::UNBOUNDED>solve_code;
  }

  bool IsProblemUnbounded() const {
    assert( IsSolStatusRetrieved() );
    return sol::INFEASIBLE<solve_code &&
        sol::UNBOUNDED>=solve_code;
  }

  bool IsSolStatusRetrieved() const {
    return sol::NOT_CHECKED!=solve_code;
  }

  struct Stats {
    steady_clock::time_point time;
    double setup_time;
    double solution_time;
  };
  Stats stats;


  /////////////////////////////// SOME MATHS ////////////////////////////////
  static bool float_equal(double a, double b) {           // ??????
    return std::fabs(a-b) < 1e-8*std::max(std::fabs(a), std::fabs(b));
  }

  bool IsFinite(double n) const {
    return n>MP_DISPATCH( MinusInfinity() ) &&
        n<MP_DISPATCH( Infinity() );
  }
  static double Infinity() { return std::numeric_limits<double>::infinity(); }
  static double MinusInfinity() { return -Infinity(); }

public:
  using Solver::add_to_long_name;
  using Solver::add_to_version;
  using Solver::set_option_header;
  using Solver::add_to_option_header;

protected:
  void HandleSolution(int status, fmt::CStringRef msg,
      const double *x, const double *y, double obj) {
    GetCQ().HandleSolution(status, msg, x, y, obj);
  }

  void HandleFeasibleSolution(fmt::CStringRef msg,
      const double *x, const double *y, double obj) {
    GetCQ().HandleFeasibleSolution(msg, x, y, obj);
  }

  ArrayRef<int> ReadSuffix(const SuffixDef<int>& suf) {
    return GetCQ().ReadSuffix(suf);
  }

  ArrayRef<double> ReadSuffix(const SuffixDef<double>& suf) {
    return GetCQ().ReadSuffix(suf);
  }

  /// Record suffix values which are written into .sol
  /// by HandleSolution()
  /// Does nothing if vector empty
  void ReportSuffix(const SuffixDef<int>& suf,
                    const std::vector<int>& values) {
    GetCQ().ReportSuffix(suf, values);
  }

  void ReportSuffix(const SuffixDef<double>& suf,
    const std::vector<double>& values) {
    GetCQ().ReportSuffix(suf, values);
  }

private:
  ///////////////////////// STORING SOLUTON STATUS //////////////////////
  int solve_code=sol::NOT_CHECKED;
  std::string solve_status;


protected:
  ///////////////////////////// OPTIONS /////////////////////////////////
  template <class Value>
  class StoredOption : public mp::TypedSolverOption<Value> {
    Value& value_;
  public:
    using value_type = Value;
    StoredOption(const char *name_list, const char *description,
                 Value& v, ValueArrayRef values = ValueArrayRef())
      : mp::TypedSolverOption<Value>(name_list, description, values), value_(v) {}

    void GetValue(Value &v) const override { v = value_; }
    void SetValue(typename internal::OptionHelper<Value>::Arg v) override
    { value_ = v; }
  };

  /// Solver options accessor, facilitates calling
  /// backend_.Get/SetSolverOption()
  template <class Value, class Index>
  class SolverOptionAccessor {
    using Backend = Impl;
    Backend& backend_;
  public:
    using value_type = Value;
    using index_type = Index;
    SolverOptionAccessor(Backend& b) : backend_(b) { }
    /// Options setup
    Value get(const SolverOption& , Index i) const {
      Value v;
      backend_.GetSolverOption(i, v);
      return v;
    }
    void set(const SolverOption& ,
             typename internal::OptionHelper<Value>::Arg v,
             Index i) {
      backend_.SetSolverOption(i, v); }
  };

  template <class ValueType, class KeyType>
  class ConcreteOptionWrapper :
      public Solver::ConcreteOptionWithInfo<
      SolverOptionAccessor<ValueType, KeyType>, ValueType, KeyType> {

    using COType = Solver::ConcreteOptionWithInfo<
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

public:
  using Solver::AddOption;
  using Solver::AddOptionSynonymsFront;
  using Solver::AddOptionSynonymsBack;
  using Solver::AddOptionSynonym_OutOfLine;
  using Solver::FindOption;


  /// Simple stored option referencing a variable
  template <class Value>
  void AddStoredOption(const char *name, const char *description,
                       Value& value, ValueArrayRef values = ValueArrayRef()) {
    AddOption(Solver::OptionPtr(
                new StoredOption<Value>(
                  name, description, value, values)));
  }

  /// Adding solver options of types int/double/string/...
  /// The type is deduced from the two last parameters min, max
  /// (currently unused otherwise - TODO)
  /// If min/max omitted, assume ValueType=std::string
  /// Assumes existence of Impl::Get/SetSolverOption(KeyType, ValueType(&))
  template <class KeyType, class ValueType=std::string>
  void AddSolverOption(const char *name_list, const char *description,
                       KeyType k,
                       ValueType , ValueType ) {
    AddOption(Solver::OptionPtr(
                new ConcreteOptionWrapper<
                ValueType, KeyType>(
                  (Impl*)this, name_list, description, k)));
  }

  template <class KeyType, class ValueType = std::string>
  void AddSolverOption(const char* name_list, const char* description,
    KeyType k) {
    AddOption(Solver::OptionPtr(
      new ConcreteOptionWrapper<
      ValueType, KeyType>(
        (Impl*)this, name_list, description, k)));
  }

  /// TODO use vmin/vmax or rely on solver raising error?
  /// TODO also with ValueTable, deduce type from it
  template <class KeyType, class ValueType = std::string>
  void AddSolverOption(const char* name_list, const char* description,
      KeyType k, ValueArrayRef values, ValueType defaultValue) {
    internal::Unused(defaultValue);
    AddOption(Solver::OptionPtr(
      new ConcreteOptionWrapper<
      ValueType, KeyType>(
        (Impl*)this, name_list, description, k, values)));
  }

  void ReplaceOptionDescription(const char* name, const char* desc) {
    auto pOption = FindOption(name);
    assert(pOption);
    pOption->set_description(desc);
  }

private:
  struct Options {
    int importPriorities_=1;
  };
  Options storedOptions_;

protected:
  void InitStandardOptions() {
    if (IMPL_HAS_STD_FEATURE( VarPriorities ))
      AddStoredOption("mip:priorities priorities",  // CP has it too
        "0/1*: Whether to read the branch and bound priorities from the"
        " .priority suffix.",
        storedOptions_.importPriorities_);
  }

  int priorities() const { return storedOptions_.importPriorities_; }

  void InitCustomOptions() { }

  ////////////////////////////////////////////////////////////
  /////////////// OPTIONAL STANDARD FEATURES /////////////////
  /////// Override in the Impl for standard operations ///////
  ////////////////////////////////////////////////////////////
protected:
  /**
  * Set branch and bound priority
  **/
  DEFINE_STD_FEATURE( VarPriorities, true ) // believe true for most
  void VarPriorities(ArrayRef<int>) {
    throw MakeUnsupportedError("BasicBackend::VarPriorities");
  }

  void ObjPriorities(ArrayRef<int>) {
    throw MakeUnsupportedError("BasicBackend::ObjPriorities");
  }
  void ObjWeights(ArrayRef<double>) {
    throw MakeUnsupportedError("BasicBackend::ObjWeights");
  }
  void ObjAbsTol(ArrayRef<double>) {
    throw MakeUnsupportedError("BasicBackend::ObjAbsTol");
  }
  void ObjRelTol(ArrayRef<double>) {
    throw MakeUnsupportedError("BasicBackend::ObjRelTol");
  }


  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////// STANDARD SUFFIXES ///////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

private:
  const SuffixDef<int> suf_varpriority = { "priority", suf::VAR | suf::INPUT };

  const SuffixDef<int> suf_objpriority = { "objpriority", suf::OBJ | suf::INPUT };
  const SuffixDef<double> suf_objweight = { "objweight", suf::OBJ | suf::INPUT };
  const SuffixDef<double> suf_objabstol = { "objabstol", suf::OBJ | suf::INPUT };
  const SuffixDef<double> suf_objreltol = { "objreltol", suf::OBJ | suf::INPUT };
};

}  // namespace mp

#endif  // BACKEND_H_
