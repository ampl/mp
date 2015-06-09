/*
 ASL expression visitor

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

#ifndef MP_ASL_ASLEXPR_VISITOR_H_
#define MP_ASL_ASLEXPR_VISITOR_H_

#include "mp/basic-expr-visitor.h"
#include "asl/aslexpr.h"

namespace mp {
namespace asl {

// ASL expression visitor.
template <typename Impl, typename Result>
class ExprVisitor :
    public BasicExprVisitor<Impl, Result, internal::ExprTypes> {};

// Returns true iff e is a zero constant.
inline bool IsZero(NumericExpr e) {
  NumericConstant c = Cast<NumericConstant>(e);
  return c && c.value() == 0;
}
}  // namespace asl
}  // namespace mp

#endif  // MP_ASL_ASLEXPR_VISITOR_H_
