/*
 Solver features

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

#ifndef TESTS_SOLVER_FEATURE_H_
#define TESTS_SOLVER_FEATURE_H_

// Solver features.
namespace feature {
enum Feature {
  FLOAT_CONST    = 0x001,
  DIV            = 0x002,
  POW            = 0x004,
  SQRT           = 0x008,
  LOG            = 0x010,
  EXP            = 0x020,
  TRIGONOMETRIC  = 0x040,
  HYPERBOLIC     = 0x080,
  PLTERM         = 0x100,
  MULTIOBJ       = 0x200,
  INITIAL_VALUES = 0x400,
  ALL            = 0xfff
};
}

#endif  // TESTS_SOLVER_FEATURE_H_
