/*
 IBM/ILOG CP direct solver interface for AMPL.

 Copyright (C) 2020 AMPL Optimization Inc

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

 Author: Gleb Belov <gleb.belov@monash.edu>
 */

#ifndef MP_INTERFACES_ILOGCP_H_
#define MP_INTERFACES_ILOGCP_H_

#include "ilogcp.h"

namespace mp {

/// TODO Currently we have to derive from IlogCPSolver
/// due to dependencies for NL file reading, interrupts, etc
/// TODO Expand for options modification, callbacks, etc
/// TODO(?) Common options, such as N threads, tolerances, etc
class IlogCPDirect : public IlogCPSolver {
public:
  /// Should be called before every model modification phase
  using IlogCPSolver::InitProblemModificationPhase;
  /// Model modification methods
  using IlogCPSolver::AddVariables;
  using IlogCPSolver::AddCommonExpressions;
  using IlogCPSolver::AddObjectives;
  using IlogCPSolver::AddAlgebraicConstraints;
  using IlogCPSolver::AddLogicalConstraints;
  /// Should be called to finalize each model modification phase
  using IlogCPSolver::FinishProblemModificationPhase;

  /// Solving
  using IlogCPSolver::Resolve;
};

}

#endif  // MP_INTERFACES_ILOGCP_H_
