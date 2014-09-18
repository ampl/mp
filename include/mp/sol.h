/*
 .sol format support.

 .sol is a format for representing solutions of mathematical optimization
 problems. It is described in the technical report "Hooking Your Solver
 to AMPL" (http://www.ampl.com/REFS/hooking2.pdf).

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

#ifndef MP_SOL_H_
#define MP_SOL_H_

#include "mp/format.h"
#include "mp/problem-base.h"

namespace mp {

struct SolutionRef {
  fmt::StringRef message;
  int num_options;
  int options[MAX_NL_OPTIONS];
  int num_vars;
  const double *values;
};

// Writes a solution in .sol format.
void WriteSol(fmt::StringRef filename, const SolutionRef &sol);
}

#endif  // MP_SOL_H_
