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

#include "mp/clock.h"
#include "mp/convert/model.h"
#include "mp/solver.h"
#include "mp/convert/constraint_keeper.h"
#include "mp/convert/std_constr.h"

namespace mp {

/// Basic backend wrapper.
/// Used by converter to directly access a solver.
/// The basic wrapper provides common functionality: option handling
/// and placeholders for solver API
template <class Impl, class Model = BasicModel< > >
class BasicBackend : public BasicConstraintAdder,
    public SolverImpl<Model>
{
public:
  BasicBackend() :
    SolverImpl<Model>(
      MP_DISPATCH( GetAMPLSolverName() ),
      MP_DISPATCH( GetAMPLSolverLongName() ),
      MP_DISPATCH( Date() ),
      MP_DISPATCH( Flags() ) )
  {
  }

  ~BasicBackend() { }

  void OpenSolver() { }
  void CloseSolver() { }
  void InitOptions() { }

  /// Default metadata
  static const char* GetSolverName() { return "SomeSolver"; }
  static std::string GetSolverVersion() { return "-285.68.53"; }
  static const char* GetAMPLSolverName() { return "solverdirect"; }
  static const char* GetAMPLSolverLongName() { return nullptr; }
  static const char* GetBackendName()    { return "BasicBackend"; }
  static const char* GetBackendLongName() { return nullptr; }
  static long Date() { return MP_DATE; }

  /// Default flags
  int Flags() {
    int flg=0;
    if (MP_DISPATCH( IfMultipleSol() ))
      flg |= Solver::MULTIPLE_SOL;
    if (MP_DISPATCH( IfMultipleObj() ))
      flg |= Solver::MULTIPLE_OBJ;
    return flg;
  }
  static bool IfMultipleSol() { return false; }
  static bool IfMultipleObj() { return false; }

  void InitializationAfterOpeningSolver() {
    MP_DISPATCH( InitNamesAndVersion() );
    MP_DISPATCH( InitOptions() );
  }
  void InitNamesAndVersion() {
    auto name = MP_DISPATCH( GetSolverName() );
    auto version = MP_DISPATCH( GetSolverVersion() );
    this->set_long_name( fmt::format("{} {}", name, version ) );
    this->set_version( fmt::format("AMPL/{} Optimizer [{}]",
                       name, version ) );
  }


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

  /// TODO why is this here?
  class LinearObjective {
    obj::Type sense_;
    std::vector<double> coefs_;
    std::vector<int> vars_;
  public:
    template <class CoefVec=std::initializer_list<double>,
              class VarVec=std::initializer_list<int> >
    LinearObjective(obj::Type s, CoefVec&& c, VarVec&& v) :
      sense_(s),
      coefs_(std::forward<CoefVec>(c)), vars_(std::forward<VarVec>(v)) { }
    obj::Type get_sense() const { return sense_; }
    int get_num_terms() const { assert(check()); return (int)vars_.size(); }
    bool check() const { return coefs_.size()==vars_.size(); }
    const std::vector<double>& get_coefs() const { return coefs_; }
    const std::vector<int>& get_vars() const { return vars_; }
  };

  void AddObjective(typename Model::Objective obj) {
    if (obj.nonlinear_expr()) {
      MP_DISPATCH( AddGeneralObjective( obj ) );
    } else {
      LinearExprUnzipper leu(obj.linear_expr());
      MP_DISPATCH( AddLinearObjective( { obj.type(),
                                         std::move(leu.c_), std::move(leu.v_) } ) );
      // TODO quadratics like in AddAlgebraicConstraint
    }
    }
  void AddLinearObjective( const LinearObjective& ) {
    throw MakeUnsupportedError("BasicBackend::AddLinearObjective");
  }
  void AddGeneralObjective(typename Model::Objective obj) {
    throw MakeUnsupportedError("BasicBackend::AddGeneralObjective");
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

  void AddGeneralConstraint(typename Model::AlgebraicCon con) {
    throw MakeUnsupportedError("BasicBackend::AddGeneralConstraint");
  }

  ////////////////// Some basic custom constraints /////////////////
  USE_BASE_CONSTRAINT_HANDLERS(BasicConstraintAdder)

  /// Optionally exclude LFDs from being posted,
  /// then all those are converted to LinearConstraint's first
  ACCEPT_CONSTRAINT(LinearDefiningConstraint, NotAccepted)
  void AddConstraint(const LinearDefiningConstraint& ldc) {
    MP_DISPATCH( AddConstraint(ldc.to_linear_constraint()) );
  }

  ACCEPT_CONSTRAINT(LinearConstraint, Recommended)
  /// TODO Attributes (lazy/user cut, etc)
  void AddConstraint(const LinearConstraint& ldc) {
    throw MakeUnsupportedError("BasicBackend::AddLinearConstraint");
  }


  void Solve(Model &p, SolutionHandler &sh) { Resolve(p, sh); }

  void Resolve(Model& p, SolutionHandler &sh) {
    MP_DISPATCH( SetInterrupter(MP_DISPATCH( interrupter() )) );

    stats.setup_time = GetTimeAndReset(stats.time);
    MP_DISPATCH( DoOptimize() );
    stats.solution_time = GetTimeAndReset(stats.time);

    // Convert solution status.
    int solve_code = 0;
    std::string status = MP_DISPATCH(
        ConvertSolutionStatus(*MP_DISPATCH( interrupter() ), solve_code) );

    fmt::MemoryWriter writer;
    writer.write("{}: {}\n", MP_DISPATCH( long_name() ), status);
    double obj_value = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> solution, dual_solution;
    if (solve_code < sol::INFEASIBLE) {
      MP_DISPATCH( PrimalSolution(solution) );

      if (MP_DISPATCH( IsMIP() )) {
        writer << MP_DISPATCH( NodeCount() ) << " nodes, ";
      } else {                                    // Also for QCP
        MP_DISPATCH( DualSolution(dual_solution) );
      }
      writer << MP_DISPATCH( Niterations() ) << " iterations";

      if (MP_DISPATCH( NumberOfObjectives() ) > 0) {
        writer.write(", objective {}",
                     MP_DISPATCH( FormatObjValue(MP_DISPATCH( ObjectiveValue() )) ));
      }

//      p.AddIntSuffix("toy_var_suffix", suf::VAR | suf::OUTPUT | suf::OUTONLY, 0);

//      /// To be separated into a separate method, like FillOutputSuffixes()
//      /// The backend writer should see a very simple interface, like
//      /// AddSuffix("name", kind, descr, table, {input_fn, output_fn});
//      MutIntSuffix toy_var_suffix =
//          Cast<MutIntSuffix>( p. suffixes(suf::VAR). Find("toy_var_suffix") );
//      assert(toy_var_suffix);
//      if (toy_var_suffix) {
//        for (int i = 0, n = MP_DISPATCH( NumberOfVariables() ); i < n; ++i) {
//          if (i % 2 == 0)
//            toy_var_suffix.set_value(i, 900+i);
//        }
//      }
    }

    sh.HandleSolution(solve_code, writer.c_str(),
        solution.empty() ? 0 : solution.data(),
        dual_solution.empty() ? 0 : dual_solution.data(), obj_value);

    double output_time = GetTimeAndReset(stats.time);

    if (MP_DISPATCH( timing() )) {
      MP_DISPATCH( Print("Setup time = {:.6f}s\n"
            "Solution time = {:.6f}s\n"
            "Output time = {:.6f}s\n",
            stats.setup_time, stats.solution_time, output_time) );
    }
  }



  /////////////////////////////// SERVICE STUFF ///////////////////////////////////
  ///
  /////////////////////////////////////////////////////////////////////////////////

  struct Stats {
    steady_clock::time_point time;
    double setup_time;
    double solution_time;
  };
  Stats stats;


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
  using Solver::add_to_option_header;

  ///////////////////////////// OPTIONS /////////////////////////////////
  /// TODOs
  /// - hide all Solver stuff behind an abstract interface
protected:

  using Solver::AddOption;

  template <class Value>
  class StoredOption : public mp::TypedSolverOption<Value> {
    Value& value_;
  public:
    using value_type = Value;
    StoredOption(const char *name, const char *description,
        Value& v, ValueArrayRef values = ValueArrayRef())
    : mp::TypedSolverOption<Value>(name, description, values), value_(v) {}

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
                          KeyType k) :
      COType(name, description, &soa_, &SOAType::get, &SOAType::set, k),
      soa_(*impl_)
    { }
  };

  public:

  /// Simple stored option referencing a variable
  template <class Value>
  void AddStoredOption(const char *name, const char *description,
                 Value& value, ValueArrayRef values = ValueArrayRef()) {
    AddOption(Solver::OptionPtr(
                      new StoredOption<Value>(
            name, description, value, values)));
  }

  /// Adding solver options of types int/double/string/...
  /// Assumes existence of Impl::Get/SetSolverOption(KeyType, ValueType(&))
  template <class KeyType, class ValueType=std::string>
  void AddSolverOption(const char *name, const char *description,
                       KeyType k,
                       /// If min/max omitted, assume ValueType=std::string
                       ValueType vMin={}, ValueType vMax={}) {
    AddOption(Solver::OptionPtr(
                      new ConcreteOptionWrapper<
                      ValueType, KeyType>(
                        (Impl*)this, name, description, k)));
  }
  /// TODO use vmin/vmax or rely on solver raising error?
  /// TODO also with ValueTable, deduce type from it


};


}  // namespace mp

#endif  // BACKEND_H_
