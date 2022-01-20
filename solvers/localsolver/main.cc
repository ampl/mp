/*
 AMPL solver interface to LocalSolver.

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

#include "localsolver/localsolver.h"
#include <lsversion.h>

extern "C" int MP_RunSolver(char **argv, const char *lic) {
  try {
// Not there in 10.5
//    if (lic)
//      localsolver::ls_set_license_content(lic);
    return mp::SolverApp<mp::LocalSolver>().Run(argv, mp::READ_BOUNDS_FIRST);
  } catch (const std::exception &e) {
    fmt::print(stderr, "Error: {}\n", e.what());
  }
  return 1;
}

#ifndef MP_NOMAIN
int main(int, char **argv) {
  return MP_RunSolver(argv, 0);
}
#endif
