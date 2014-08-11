/*
 An AMPL function library that provides functions used by
 constraint programming solvers.

 Copyright (C) 2013 AMPL Optimization Inc

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

#include <cstring>
#include <limits>
#include "mp/format.h"
#include "funcadd.h"

namespace {
void error(arglist *al, fmt::StringRef message) {
  al->Errmsg = static_cast<char*>(al->AE->Tempmem(al->TMI, message.size() + 1));
  std::strcpy(al->Errmsg, message.c_str());
}

double element(arglist *al) {
  if (al->derivs) {
    for (int i = 0, n = al->n; i < n; ++i)
      al->derivs[i] = std::numeric_limits<double>::quiet_NaN();
    if (al->hes) {
      error(al, "derivatives are not provided");
      return 0;
    }
  }
  double index = al->ra[al->n - 1];
  int int_index = static_cast<int>(index);
  if (int_index != index || int_index < 0 || int_index >= al->n - 1) {
    fmt::Writer message;
    message.write("invalid index {}", index);
    error(al, message.c_str());
    return 0;
  }
  return al->ra[int_index];
}

double in_relation(arglist *al) {
  static const char error[] = "can't evaluate in_relation";
  al->Errmsg = static_cast<char*>(al->AE->Tempmem(al->TMI, sizeof(error)));
  std::strcpy(al->Errmsg, error);
  return 0;
}
}

extern "C" void funcadd_ASL(AmplExports *ae) {
  enum {MIN_ELEMENT_ARGS = 2};
  ae->Addfunc("element", element,
      FUNCADD_REAL_VALUED, -MIN_ELEMENT_ARGS - 1, 0, ae);

  // in_relation takes at least 1 argument, a variable.
  enum {MIN_IN_RELATION_ARGS = 1};
  ae->Addfunc("in_relation", in_relation,
      FUNCADD_REAL_VALUED, -MIN_IN_RELATION_ARGS - 1, 0, ae);
}
