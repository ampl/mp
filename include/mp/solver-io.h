#ifndef SOLVERIO_H
#define SOLVERIO_H

/**
 * Support for .nl/.sol/console I/O requiring a (Basic)Solver
 */

#include <functional>

#include "mp/arrayref.h"

#include "mp/solver-base.h"

#include "mp/nl-reader.h"
#include "mp/sol.h"

namespace mp {

/// Adapts a solution for WriteSol.
template <typename ProblemBuilder>
class SolutionAdapter {
 private:
  int status_;
  ProblemBuilder *builder_;
  const char *message_;
  mp::ArrayRef<int> options_;
  mp::ArrayRef<double> values_;
  mp::ArrayRef<double> dual_values_;
  int objno_;

 public:
  SolutionAdapter(int status, ProblemBuilder *pb, const char *message,
                  mp::ArrayRef<int> options, mp::ArrayRef<double> values,
                  mp::ArrayRef<double> dual_values,
                  int on)
    : status_(status), builder_(pb), message_(message), options_(options),
      values_(values), dual_values_(dual_values), objno_(on) {}

  int status() const { return status_; }

  const char *message() const { return message_; }

  int num_options() const { return static_cast<int>(options_.size()); }
  int option(int index) const { return options_[index]; }

  int num_values() const { return static_cast<int>(values_.size()); }
  double value(int index) const { return values_[index]; }

  int num_dual_values() const { return static_cast<int>(dual_values_.size()); }
  double dual_value(int index) const { return dual_values_[index]; }

  int objno() const { return objno_; }

  int num_vars() const { return builder_->num_vars(); }
  int num_algebraic_cons() const { return builder_->num_algebraic_cons(); }

  const typename ProblemBuilder::SuffixSet *suffixes(suf::Kind kind) const {
    return builder_ ? &builder_->suffixes(kind) : 0;
  }
};

/// The default .sol file writer.
class SolFileWriter {
 public:
  template <typename Solution>
  void Write(fmt::CStringRef filename, const Solution &sol) {
    WriteSolFile(filename, sol);
  }
};

/// A solution writer.
/// @param Solver: optimization solver class
/// @param ProblemBuilder: suffix handler
/// @param Writer: .sol writer
template <typename Solver,
          typename ProblemBuilder,
          typename Writer = SolFileWriter>
class SolutionWriterImpl :
    private Writer, public SolutionHandler {
 private:
  std::string stub_;
  std::string overrideStub_;
  Solver &solver_;
  ProblemBuilder &builder_;

  ArrayRef<int> options_;

  /// The number of feasible solutions found.
  int num_solutions_;

 protected:
  Solver &solver() { return solver_; }
  ProblemBuilder &builder() { return builder_; }
  const std::string &stub() const { return stub_; }

 public:
  
  SolutionWriterImpl(fmt::StringRef stub, Solver &s, ProblemBuilder &b,
                 ArrayRef<int> options = mp::ArrayRef<int>(0, 0))
    : stub_(stub.to_string()), overrideStub_(), solver_(s), builder_(b),
      options_(options), num_solutions_(0) {}
  /// Returns the .sol writer.
  Writer &sol_writer() { return *this; }

  /// Writes an intermediate solution to a .sol file.
  void HandleFeasibleSolution(int status, fmt::CStringRef message,
        const double *values, const double *dual_values, double);

  /// Deprecated: no status
  void HandleFeasibleSolution(fmt::CStringRef message,
        const double *values, const double *dual_values, double obj) {
    HandleFeasibleSolution(sol::UNCERTAIN, message,
                           values, dual_values, obj);
  }

  /// Writes the solution to a .sol file.
  void HandleSolution(int status, fmt::CStringRef message,
        const double *values, const double *dual_values, double);

  virtual void OverrideSolutionFileName(const std::string& fileName) {
    overrideStub_ = fileName;
  };
};

/// Convenience typedef for APIs using SolverImpl<> or similar
template <class Solver, class Writer = SolFileWriter>
using SolutionWriter = SolutionWriterImpl<
                           Solver,
                           typename Solver::ProblemBuilder,
                           Writer>;

template <typename Solver, typename PB, typename Writer>
void SolutionWriterImpl<Solver, PB, Writer>::
HandleFeasibleSolution(
    int status, fmt::CStringRef message, const double *values,
    const double *dual_values, double) {
  ++num_solutions_;
  const char *solution_stub = solver_.solution_stub();
  if (!*solution_stub)
    return;
  SolutionAdapter<PB> sol(
        status, &builder_, message.c_str(), options_,
        MakeArrayRef(values, values ? builder_.num_vars() : 0),
        MakeArrayRef(dual_values,
                     dual_values ? builder_.num_algebraic_cons() : 0),
        solver_.objno_used());
  fmt::MemoryWriter filename;
  filename << solution_stub << num_solutions_ << ".sol";
  this->Write(filename.c_str(), sol);
}

template <typename Solver, typename PB, typename Writer>
void SolutionWriterImpl<Solver, PB, Writer>::HandleSolution(
    int status, fmt::CStringRef message, const double *values,
    const double *dual_values, double) {
  if (solver_.need_multiple_solutions()) {
    auto kindP = mp::suf::Kind( suf::PROBLEM | suf::OUTPUT | suf::OUTONLY );
    auto kindO = mp::suf::Kind( suf::OBJ | suf::OUTPUT | suf::OUTONLY );
    builder_.AddIntSuffix("nsol", kindP, 0).
        SetValue(0, num_solutions_);
    builder_.AddIntSuffix("nsol", kindO, 0).
        SetValue(0, num_solutions_);
    builder_.AddIntSuffix("npool", kindP, 0).
        SetValue(0, num_solutions_);
    builder_.AddIntSuffix("npool", kindO, 0).
        SetValue(0, num_solutions_);
  }
  SolutionAdapter<PB> sol(
        status, &builder_, message.c_str(), options_,
        MakeArrayRef(values, values ? builder_.num_vars() : 0),
        MakeArrayRef(dual_values,
                     dual_values ? builder_.num_algebraic_cons() : 0),
        solver_.objno_used());

  std::string solFilePath;
  if (!overrideStub_.empty()) { 
    // Check if its absolute or relative
    if ((overrideStub_.length() >= 2) &&// cannot be absolute otherwise 
       ((overrideStub_[0] == '/') || // linux abs (/ ...)
        (overrideStub_[1] == ':') ))  {// windows abs (x: ///)
        solFilePath = overrideStub_; // if absolute
    }
    else { // if relative
      size_t pos = stub_.find_last_of("\\/"); //.find parent dir
      auto dir = (std::string::npos == pos) ? "" : stub_.substr(0, pos+1);
      solFilePath = dir + overrideStub_;
    }
  }
  else { solFilePath = stub_ + ".sol"; }
  
  this->Write(solFilePath, sol);
}


namespace internal {

/// An .nl handler for an application using BasicSolver.
/// The full definition separates parameters
/// Solver, ProblemBuilder, and NLProblemBuilder.
template <class Solver,
          class ProblemBuilder,
          class NLProblemBuilder>
class SolverNLHandlerImpl : public NLProblemBuilder {
private:
  Solver &solver_;
  int num_options_;
  int options_[MAX_AMPL_OPTIONS];
  std::function<void()> after_header_;

  typedef NLProblemBuilder Base;

public:
  SolverNLHandlerImpl(ProblemBuilder &pb,
                      Solver &s,
                      std::function<void()> after_h = {})
    : Base(pb), solver_(s), num_options_(0),
      after_header_(after_h)
  { }

  int num_options() const { return num_options_; }
  const int *options() const { return options_; }

  int objno() const override { return solver_.objno_specified(); }
  bool multiobj() const override { return solver_.multiobj(); }
  void notify_obj_added() const override { solver_.notify_obj_added(); }

  void OnHeader(const NLHeader &h);
};


/// Perform demo version checks if necessary.
void CheckDemoVersion(const NLHeader &h);


/// A shorthand typedef for APIs using SolverImpl or similar
template <class Solver>
using SolverNLHandler = SolverNLHandlerImpl<
                              Solver,
                              typename Solver::ProblemBuilder,
                              typename Solver::NLProblemBuilder>;

template <typename Solver, typename PB, typename NLPB>
void SolverNLHandlerImpl<Solver, PB, NLPB>::OnHeader(const NLHeader &h) {
  num_options_ = h.num_ampl_options;
  std::copy(h.ampl_options, h.ampl_options + num_options_, options_);
  if (after_header_) {
    solver_.notify_start_opts();
    after_header_();
  }
  solver_.notify_end_opts();
  /// Clarify objectives
  int objno = solver_.objno_specified();
  if (objno > h.num_objs && solver_.is_objno_specified())
    throw InvalidOptionValue("objno", objno,
                             fmt::format("expected value between 0 and {}", h.num_objs));
  Base::OnHeader(h);
#ifndef MP_DATE
  CheckDemoVersion(h);
#endif
}


/// Prints a solution to stdout.
void PrintSolution(const double *values, int num_values, const char *name_col,
                   const char *value_col, NameProvider &np);

/// Solution handler for a solver application.
/// Extends SolutionWriterImpl by consideration of
/// the -AMPL switch and the \a wantsol option.
template <class Solver,
          class ProblemBuilder,
          class Writer = SolFileWriter>
class AppSolutionHandlerImpl :
    public SolutionWriterImpl<Solver, ProblemBuilder, Writer> {
 private:
  unsigned banner_size_;

 public:
  AppSolutionHandlerImpl(fmt::StringRef stub, Solver &s,
                    ProblemBuilder &b,
                    ArrayRef<int> options, unsigned banner_size)
  : SolutionWriterImpl<Solver, ProblemBuilder, Writer>
      (stub, s, b, options),  banner_size_(banner_size) {}

  /// Specialize HandleSolution()
  void HandleSolution(int status, fmt::CStringRef message, const double *values,
                      const double *dual_values, double obj_value);
};


/// Convenience typedef for APIs using SolverImpl<> or similar
template <class Solver, class Writer = SolFileWriter>
using AppSolutionHandler = AppSolutionHandlerImpl<
                             Solver,
                             typename Solver::ProblemBuilder,
                             Writer>;


template <typename Solver, typename PB, typename Writer>
void AppSolutionHandlerImpl<Solver, PB, Writer>::HandleSolution(
    int status, fmt::CStringRef message, const double *values,
    const double *dual_values, double obj_value) {
  Solver &solver = this->solver();
  int wantsol = solver.wantsol();
  if (solver.ampl_flag() || (wantsol & Solver::WRITE_SOL_FILE) != 0) {
    // "Erase" the banner so that it is not duplicated when printing
    // the solver message.
    if (solver.ampl_flag() && banner_size_ != 0) {
      fmt::MemoryWriter w;
      w << fmt::pad("", banner_size_, '\b');
      solver.Print("{}", w.c_str());
    }
    SolutionWriterImpl<Solver, PB, Writer>::HandleSolution(
          status, message, values, dual_values, obj_value);
  }
  if (solver.ampl_flag())
    return;
  if ((wantsol & Solver::SUPPRESS_SOLVER_MSG) == 0)
    solver.Print("{}\n", message.c_str() + banner_size_);
  using internal::PrintSolution;
  if ((wantsol & Solver::PRINT_SOLUTION) != 0) {
    int num_vars = this->builder().num_vars();
    NameProvider np(this->stub() + ".col", "_svar", num_vars);
    PrintSolution(values, num_vars, "variable", "value", np);
  }
  if ((wantsol & Solver::PRINT_DUAL_SOLUTION) != 0) {
    int num_cons = this->builder().num_algebraic_cons();
    NameProvider np(this->stub() + ".row", "_scon", num_cons);
    PrintSolution(dual_values, num_cons, "constraint", "dual value", np);
  }
}

}  // namespace internal

} // namespace mp

#endif // SOLVERIO_H
