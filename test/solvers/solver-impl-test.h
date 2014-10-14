/*
 Solver implementation test

 Copyright (C) 2012 AMPL Optimization Inc

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

 Author: Victor Zverovich
 */

#ifndef TESTS_SOLVER_TEST_H_
#define TESTS_SOLVER_TEST_H_

#include "asl/aslbuilder.h"
#include "asl/aslsolver.h"
#include "../solution-handler.h"
#include "gtest/gtest.h"

#ifdef HAVE_THREADS
# include <thread>
#endif

// Solver features.
namespace feature {
enum Feature {
  FLOAT_CONST = 0x01,
  DIV         = 0x02,
  POW         = 0x04,
  HYPERBOLIC  = 0x08,
  SQRT        = 0x10,
  LOG         = 0x20,
  EXP         = 0x40,
  PLTERM      = 0x80,
  ALL         = 0xff
};
}

#ifdef MP_USE_UNIQUE_PTR
typedef std::unique_ptr<mp::ASLSolver> SolverPtr;
#else
typedef std::auto_ptr<mp::ASLSolver> SolverPtr;
#endif

struct SolverTestParam {
  typedef SolverPtr (*CreateSolver)();

  CreateSolver create_solver;
  unsigned features;

  SolverTestParam(CreateSolver cs, unsigned features)
  : create_solver(cs), features(features) {}
};

// Define a dummy operator<< because otherwise GTest produces uninitialized
// memory access under valgrind when trying to format SolverTestParam.
inline std::ostream &operator<<(std::ostream &os, const SolverTestParam &) {
  return os << "[SolverTestParam]";
}

// Abstract solver test.
class SolverImplTest
    : public ::testing::TestWithParam<SolverTestParam>,
      public mp::internal::ASLBuilder {
 protected:
  SolverPtr solver_;
  unsigned features_;
  mp::Variable x;
  mp::Variable y;
  mp::Variable z;

  FMT_DISALLOW_COPY_AND_ASSIGN(SolverImplTest);

  mp::NumericConstant MakeConst(double value) {
    return MakeNumericConstant(value);
  }

  bool HasFeature(feature::Feature f) const {
    return (features_ & f) != 0;
  }

  static double TestFunc(arglist *) { return 0; }

  class EvalResult {
   private:
    bool has_value_;
    double value_;
    double obj_value_;
    int solve_code_;

   public:
    explicit EvalResult(int solve_code = mp::sol::UNKNOWN)
    : has_value_(false), value_(), obj_value_(), solve_code_(solve_code) {}

    EvalResult(double value, double obj_value)
    : has_value_(true), value_(value), obj_value_(obj_value),
      solve_code_(mp::sol::UNKNOWN) {}

    bool has_value() const { return has_value_; }

    friend bool operator==(double lhs, const EvalResult &rhs) {
      if (!rhs.has_value_)
        throw std::runtime_error("no value");
      return lhs == rhs.value_;
    }

    double obj_value() const {
      if (!has_value_)
        throw std::runtime_error("no value");
      return obj_value_;
    }
    int solve_code() const { return solve_code_; }
    void set_solve_code(int code) { solve_code_ = code; }
  };

  EvalResult Solve(mp::Problem &p);

  // Solves a problem containing a single constraint with the given
  // expression and returns the value of the variable with index 0.
  EvalResult Solve(mp::LogicalExpr e,
      int var1, int var2, int var3 = 0, bool need_result = false);

  // Evaluates a numeric expression by constructing and solving a problem.
  EvalResult Eval(mp::NumericExpr e,
      int var1 = 0, int var2 = 0, int var3 = 0) {
    return Solve(MakeRelational(mp::expr::EQ, MakeVariable(0), e),
                 var1, var2, var3, true);
  }

  // Evaluates a logical expression by constructing and solving a problem.
  EvalResult Eval(mp::LogicalExpr e,
      int var1 = 0, int var2 = 0, int var3 = 0) {
    return Eval(MakeIf(e, MakeConst(1), MakeConst(0)), var1, var2, var3);
  }

  SolveResult Solve(mp::Problem &p, const char *stub, const char *opt = 0) {
    return Solve(*solver_, p, stub, opt);
  }

  SolveResult Solve(const char *stub, const char *opt = 0) {
    mp::Problem p;
    return Solve(*solver_, p, stub, opt);
  }

 public:
  SolverImplTest();

  static SolveResult Solve(mp::ASLSolver &s,
      mp::Problem &p, const char *stub, const char *opt = 0);
};

void Interrupt();

#ifndef MP_TEST_DATA_DIR
# define MP_TEST_DATA_DIR "../data"
#endif

#endif  // TESTS_SOLVER_TEST_H_
