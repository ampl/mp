#ifndef TESTS_SOLUTION_HANDLER_H_
#define TESTS_SOLUTION_HANDLER_H_

#include <limits>

#include "mp/solver.h"

class TestSolutionHandler : public mp::BasicSolutionHandler {
 private:
  int num_vars_;
  int status_;
  std::string message_;
  std::vector<double> primal_;
  const double *dual_;
  double obj_value_;

 public:
  explicit TestSolutionHandler(int num_vars = 0)
  : num_vars_(num_vars), status_(mp::sol::UNKNOWN), primal_(0), dual_(0),
    obj_value_(std::numeric_limits<double>::quiet_NaN()) {}
  virtual ~TestSolutionHandler() {}

  int status() const { return status_; }
  double obj_value() const { return obj_value_; }
  const std::string &message() const { return message_; }
  const double *primal() const { return primal_.empty() ? 0 : &primal_[0]; }
  const double *dual() const { return dual_; }

  void HandleSolution(int status, fmt::StringRef message,
        const double *primal, const double *dual, double obj_value) {
    status_ = status;
    message_ = message;
    if (primal)
      primal_.assign(primal, primal + num_vars_);
    dual_ = dual;
    obj_value_ = obj_value;
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
