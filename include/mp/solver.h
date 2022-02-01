/*
 A mathematical optimization solver.

 Copyright (C) 2012-2021 AMPL Optimization Inc

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

#ifndef MP_SOLVER_H_
#define MP_SOLVER_H_

#include <memory>
#include <string>
#include <vector>

#include "mp/problem-builder.h"
#include "mp/solver-base.h"
#include "mp/solver-io.h"
#include "mp/solver-app.h"
#include "mp/suffix.h"

namespace mp {

/** 
  DEPRECATED.

  An old base class for solver drivers.
  */
class Solver : public BasicSolver {

public:
  /// Constructs a Solver object.
  /// date:  The solver date in YYYYMMDD format.
  /// flags: Bitwise OR of zero or more of the following values
  ///          MULTIPLE_SOL
  ///          MULTIPLE_OBJ
  Solver(fmt::CStringRef name, fmt::CStringRef long_name,
              long date, int flags) :
    BasicSolver(name, long_name, date, flags) { }


  class SuffixInfo {
   private:
    const char *name_;
    const char *table_;
    int kind_;
    int nextra_;

   public:
    SuffixInfo(const char *name, const char *table, int kind, int nextra)
      : name_(name), table_(table), kind_(kind), nextra_(nextra) {}

    const char *name() const { return name_; }
    const char *table() const { return table_; }
    int kind() const { return kind_; }
    int nextra() const { return nextra_; }
  };

  typedef std::vector<SuffixInfo> SuffixList;

private:
  SuffixList suffixes_;

public:
  /// Adds a suffix.
  void AddSuffix(const char *name, const char *table,
                 int kind, int nextra = 0) {
    suffixes_.push_back(SuffixInfo(name, table, kind, nextra));
  }

  const SuffixList &suffixes() const { return suffixes_; }

};


/// DEPRECATED.
/// Convenience template used by some APIs
/// to provide [NL]ProblemBuilderType
template <typename ProblemBuilderT>
class SolverImpl : public Solver {
 public:
  typedef ProblemBuilderT ProblemBuilder;
  typedef internal::NLProblemBuilder<ProblemBuilder> NLProblemBuilder;

  SolverImpl(fmt::CStringRef name, fmt::CStringRef long_name = 0,
             long date = 0, int flags = 0)
    : Solver(name, long_name, date, flags) {}
};


#ifdef MP_USE_UNIQUE_PTR
typedef std::unique_ptr<Solver> SolverPtr;
#else
typedef std::auto_ptr<Solver> SolverPtr;
inline SolverPtr move(SolverPtr p) { return p; }
#endif

}  // namespace mp

#endif  // MP_SOLVER_H_
