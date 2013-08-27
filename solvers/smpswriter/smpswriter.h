/*
 SMPS writer implemented as an AMPL solver.

 Copyright (C) 2013 AMPL Optimization LLC

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

#ifndef AMPL_SOLVERS_SMPSWRITER_H
#define AMPL_SOLVERS_SMPSWRITER_H

#include "solvers/util/solver.h"

namespace ampl {

class SMPSWriter : public Solver<SMPSWriter> {
 private:
  std::string ExtractScenario(std::string &name) const;

 public:
  SMPSWriter();

  void Solve(Problem &p);
};
}

#endif // AMPL_SOLVERS_SMPSWRITER_H
