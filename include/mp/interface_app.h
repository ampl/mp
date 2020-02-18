/*
 A mathematical optimization app using abstract interfaces.

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

#ifndef INTERFACE_APP_H_
#define INTERFACE_APP_H_


#include "mp/solver.h"
#include "mp/backend.h"


namespace mp {

template <class Interface>
class InterfaceApp : protected SolverApp<Interface> {
public:
  // Runs the application.
  // It processes command-line arguments and, if the NL file name (stub) is
  // specified, reads the problem, parses solver options
  // from argv and environment variables, converts and solves the problem
  // and writes solution(s).
  // argv: an array of command-line arguments terminated by a null pointer
  int RunFromNLFile(char **argv, int nl_reader_flags = 0) {
    if (!Init( argv, nl_reader_flags ))
      return GetResultCode();
    ReadNL( nl_reader_flags );
    ConvertModel();
    Resolve();
    return 0;
  }
protected:
  using TBaseClass = SolverApp<Interface>;

  using TBaseClass::Init;
  using TBaseClass::ReadNL;
  void ConvertModel() {
    ModelToBackendFeeder<Problem, Interface>
        feeder(GetProblemBuilder().problem(), GetSolver());
    feeder.PushWholeProblem();
  }
  using TBaseClass::Resolve;

  using TBaseClass::GetResultCode;
  using TBaseClass::GetProblemBuilder;
  using TBaseClass::GetSolver;
};

}  // namespace mp

#endif  // INTERFACE_APP_H_
