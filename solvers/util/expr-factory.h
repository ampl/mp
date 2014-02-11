/*
 An AMPL expression factory.

 Copyright (C) 2014 AMPL Optimization Inc

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

#ifndef SOLVERS_UTIL_EXPR_FACTORY_H_
#define SOLVERS_UTIL_EXPR_FACTORY_H_

#include "solvers/util/expr.h"
#include "solvers/util/noncopyable.h"

namespace ampl {

class NLHeader;
class ExprFactory;

namespace internal {
// Initialize ASL from an NLHeader object.
void InitASL(ASL &asl, const char *stub, const NLHeader &h);

const ASL &GetASL(const ExprFactory &ef);
}

class ExprFactory : Noncopyable {
 private:
  ASL *asl_;
  efunc *r_ops_[N_OPS];

  friend const ASL &internal::GetASL(const ExprFactory &ef) { return *ef.asl_; }

 public:
  ExprFactory(const NLHeader &h, const char *stub);
  ~ExprFactory();

  NumericConstant CreateNumericConstant(double value);
  Variable CreateVariable(int var_index);
};
}

#endif  // SOLVERS_UTIL_EXPR_FACTORY_H_
