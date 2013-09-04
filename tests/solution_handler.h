#ifndef TESTS_SOLUTION_HANDLER_H_
#define TESTS_SOLUTION_HANDLER_H_

#include <limits>

#include "solvers/util/solver.h"

class TestSolutionHandler : public ampl::SolutionHandler {
 private:
  ampl::Problem *problem_;
  std::string message_;
  double obj_value_;
  int solve_code_;
  const double *primal_;
  const double *dual_;

 public:
  TestSolutionHandler()
  : problem_(0), obj_value_(std::numeric_limits<double>::quiet_NaN()),
    solve_code_(0), primal_(0), dual_(0) {}
  virtual ~TestSolutionHandler() {}

  ampl::Problem *problem() const { return problem_; }
  int solve_code() const { return solve_code_; }
  double obj_value() const { return obj_value_; }
  const std::string &message() const { return message_; }
  const double *primal() const { return primal_; }
  const double *dual() const { return dual_; }

  void HandleSolution(ampl::Problem &p, fmt::StringRef message,
        const double *primal, const double *dual, double obj_value) {
    problem_ = &p;
    solve_code_ = p.solve_code();
    message_ = message;
    obj_value_ = obj_value;
    primal_ = primal;
    dual_ = dual;
  }
};

struct SolveResult {
  bool solved;
  double obj;
  std::string message;
  SolveResult(bool solved, double obj, const std::string &message)
  : solved(solved), obj(obj), message(message) {}
};

#endif  // TESTS_SOLUTION_HANDLER_H_
