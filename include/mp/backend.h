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
template <class Impl, class Model = BasicModel<std::allocator<char>>>
class BasicBackend : public BasicConstraintAdder,
    public SolverImpl<Model>
{
public:
  BasicBackend(fmt::CStringRef name, fmt::CStringRef longname=0,
               long date=0, int flags=0) :
    SolverImpl<Model>(name, longname, date, flags) { }

  using Variable = Problem::Variable;

  static const char* GetBackendName() { return "BasicBackend"; }

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

  class LinearObjective {
    obj::Type sense_;
    std::vector<double> coefs_;
    std::vector<int> vars_;
  public:
    template <class CoefVec=std::initializer_list<double>,
              class VarVec=std::initializer_list<int> >
    LinearObjective(obj::Type s, CoefVec&& c, VarVec&& v) :
      sense_(s), coefs_(std::move(c)), vars_(std::move(v)) { }
    obj::Type get_sense() const { return sense_; }
    int get_num_terms() const { assert(check()); return (int)vars_.size(); }
    bool check() const { return coefs_.size()==vars_.size(); }
    const std::vector<double>& get_coefs() const { return coefs_; }
    const std::vector<int>& get_vars() const { return vars_; }
  };

  void AddObjective(Problem::Objective obj) {
    if (obj.nonlinear_expr()) {
      MP_DISPATCH( AddGeneralObjective( obj ) );
    } else {
      LinearExprUnzipper leu(obj.linear_expr());
      MP_DISPATCH( AddLinearObjective( { obj.type(),
                                         std::move(leu.c_), std::move(leu.v_) } ) );
    }
    }
  void AddLinearObjective( const LinearObjective& ) {
    throw MakeUnsupportedError("BasicBackend::AddLinearObjective");
  }
  void AddGeneralObjective(Problem::Objective obj) {
    throw MakeUnsupportedError("BasicBackend::AddGeneralObjective");
  }

  void AddAlgebraicConstraint(Problem::AlgebraicCon con) {
    if (con.nonlinear_expr()) {
      MP_DISPATCH( AddGeneralConstraint( con ) );
    } else {
      LinearExprUnzipper leu(con.linear_expr());
      MP_DISPATCH( AddLinearConstraint(leu.c_.size(), leu.c_.data(), leu.v_.data(),
                                       con.lb(), con.ub()) );
    }
    }

  /// TODO Do we need ability to add several at once?
  /// TODO Attributes (lazy/user cut, etc)
  void AddLinearConstraint(int nnz, const double* c, const int* v,
                           double lb, double ub) {
    throw MakeUnsupportedError("BasicBackend::AddLinearConstraint");
  }
  void AddGeneralConstraint(Problem::AlgebraicCon con) {
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
  void AddConstraint(const LinearConstraint& ldc) {     // TODO make this form primary
    MP_DISPATCH( AddLinearConstraint(ldc.nnz(), ldc.coefs(), ldc.vars(),
                                     ldc.lb(), ldc.ub()) );
  }


  void Solve(Problem &p, SolutionHandler &sh) { Resolve(p, sh); }

  void Resolve(Problem& p, SolutionHandler &sh) {
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

  ///////////////////////////// OPTIONS /////////////////////////////////
  ///
public:
  using Solver::add_to_long_name;
  using Solver::add_to_version;
  using Solver::add_to_option_header;

protected:

  using Solver::AddOption;

  /// Simple stored option referencing a variable
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
public:
  template <class Value>
  void AddOption(const char *name, const char *description,
                 Value& value, ValueArrayRef values = ValueArrayRef()) {
    AddOption(Solver::OptionPtr(
                      new StoredOption<Value>(
            name, description, value, values)));
  }

protected:
  /// Options stored in an 'options manager'
  template <class OptionsManager>
  void AddOption(const char *name, const char *description,
                 OptionsManager& om, typename OptionsManager::index_type i,
                 ValueArrayRef values = ValueArrayRef()) {
    AddOption(Solver::OptionPtr(
                      new Solver::ConcreteOptionWithInfo<OptionsManager,
                      typename OptionsManager::value_type, typename OptionsManager::index_type>(
            name, description, &om, &OptionsManager::get, &OptionsManager::set, i, values)));
  }


  template <class Value, class Index, Index N>
  class OptionArrayManager {
    std::array<Value, N> values_;
  public:
    using value_type = Value;
    using index_type = Index;
    /// Options setup
    Value get(const SolverOption& , Index i) const { return values_.at(i); }
    void set(const SolverOption& ,
             typename internal::OptionHelper<Value>::Arg v,
             Index i) { values_.at(i) = v; }
    /// Normal getter
    const value_type& get(Index i) const { return values_.at(i); }
  };
  template <class Index, Index N>
  class OptionArrayManager<std::string, Index, N> {
    std::array<std::string, N> values_;
  public:
    using value_type = std::string;
    using index_type = Index;
    /// Options setup
    value_type get(const SolverOption& , Index i) const { return values_.at(i); }
    void set(const SolverOption& ,
             typename internal::OptionHelper<std::string>::Arg v,
             Index i) { values_.at(i) = v.data(); }
    /// Normal getter
    const value_type& get(Index i) const { return values_.at(i); }
  };

  /// Solver options
  template <class Backend, class Value, class Index>
  class SolverOptionAccessor {
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

};


}  // namespace mp

#endif  // BACKEND_H_
