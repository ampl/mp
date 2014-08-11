/*
 Second-order stochastic dominance functions for AMPL.

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

#include "ssdsolver/ssdsolver.h"
#include "funcadd.h"

#include <cstring>

#define STRINGIFY(x) #x

namespace {
const char *ssd_version(arglist *) {
  return STRINGIFY(SSDSOLVER_VERSION);
}

double ssd_uniform(arglist *al) {
  static const char error[] = "can't evaluate ssd_uniform";
  al->Errmsg = static_cast<char*>(al->AE->Tempmem(al->TMI, sizeof(error)));
  std::strcpy(al->Errmsg, error);
  return 0;
}
}

extern "C" void funcadd_ASL(AmplExports *ae) {
  ae->Addfunc("ssd_version", reinterpret_cast<rfunc>(ssd_version),
      FUNCADD_STRING_VALUED, 0, 0, ae);
  ae->Addfunc("ssd_uniform", ssd_uniform, FUNCADD_REAL_VALUED, 2, 0, ae);
}
