/*
 A program for generating AMPL declarations for the functions provided
 by the amplgsl library.

 Copyright (C) 2012 AMPL Optimization LLC

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization LLC disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#include <stdio.h>
#include "funcadd.h"

#undef printf

#define UNUSED(x) (void)(x)

/* See AddFunc in funcadd.h */
static void declare_func(const char *name, rfunc f,
    int type, int nargs, void *funcinfo, AmplExports *ae) {
  UNUSED(f);
  UNUSED(type);
  UNUSED(nargs);
  UNUSED(funcinfo);
  UNUSED(ae);
  printf("function %s%s;\n", name,
      (type & FUNCADD_RANDOM_VALUED) != 0 ? " random" : "");
}

static void dummy_at_reset(AmplExports *ae, Exitfunc *f, void *data) {
  UNUSED(ae);
  UNUSED(f);
  UNUSED(data);
}

int main() {
  AmplExports ae = {0};
  ae.Addfunc = declare_func;
  ae.AtReset = dummy_at_reset;
  printf(
      "# Automatically generated AMPL declarations for the GSL functions.\n"
      "load amplgsl.dll;\n");
  funcadd_ASL(&ae);
  return 0;
}
