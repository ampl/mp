/*
 * A program for generating AMPL declarations for the functions provided by the
 * amplgsl library.
 */

#include <stdio.h>
#include "solvers/funcadd.h"

#undef printf

#define UNUSED(x) (void)(x)

/* See AddFunc in funcadd.h */
static void declare_func(const char *name, rfunc f,
    int type, int nargs, void *funcinfo, AmplExports *ae) {
  UNUSED(f); UNUSED(type); UNUSED(nargs); UNUSED(funcinfo); UNUSED(ae);
  printf("function %s;\n", name);
}

int main() {
  AmplExports ae = {0};
  ae.Addfunc = declare_func;
  printf(
      "# Automatically generated AMPL declarations for the GSL functions.\n"
      "load libamplgsl.so;\n");
  funcadd_ASL(&ae);
  return 0;
}
