/*
 A test AMPL function library.

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

#include "solvers/funcadd.h"

#define UNUSED(x) (void)(x)

static double TestFunc(arglist *al) {
  UNUSED(al);
  return 42;
}

static int TestHandler(AmplExports *ae, TableInfo *ti) {
  UNUSED(ae);
  UNUSED(ti);
  return DB_Done;
}

extern "C" void funcadd_ASL(AmplExports *ae) {
  addfunc("foo", TestFunc, FUNCADD_REAL_VALUED, 0, 0);
  addfunc("bar", TestFunc, FUNCADD_REAL_VALUED, 0, 0);
  addfunc("bar", TestFunc, FUNCADD_REAL_VALUED, 0, 0);
  add_table_handler(TestHandler, TestHandler,
      const_cast<char*>("testhandler\n"), 0, 0);
}
