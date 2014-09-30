#ifndef TESTS_SOLUTION_HANDLER_H_
#define TESTS_SOLUTION_HANDLER_H_

#include <limits>

#include "mp/solver.h"
#include "asl/problem.h"

class TestSolutionHandler : public mp::BasicSolutionHandler {
 private:
  mp::Problem &problem_;
  std::string message_;
  double obj_value_;
  const double *primal_;
  const double *dual_;

 public:
  explicit TestSolutionHandler(mp::Problem &p)
  : problem_(p), obj_value_(std::numeric_limits<double>::quiet_NaN()),
    primal_(0), dual_(0) {}
  virtual ~TestSolutionHandler() {}

  double obj_value() const { return obj_value_; }
  const std::string &message() const { return message_; }
  const double *primal() const { return primal_; }
  const double *dual() const { return dual_; }

  void HandleSolution(fmt::StringRef message,
        const double *primal, const double *dual, double obj_value) {
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
