/*
 AMPL solver interface to IBM/ILOG CP solver via NL model file.

 Copyright (C) 2012 - 2020 AMPL Optimization Inc

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

#include "ilogcp.h"
#include "mp/ampls-c-api.h" // for CCallbacks


extern "C" int main1(int, char **argv) {
  // Solver should be destroyed after any IloException is handled.
  try {
    typedef mp::SolverApp<mp::IlogCPSolver> IlogCPApp;
    std::unique_ptr<IlogCPApp> s;
    try {
      s.reset(new IlogCPApp());
      return s->Run(argv);
    } catch (const IloException &e) {
      fmt::print(stderr, "Error: {}\n", e.getMessage());
    }
  } catch (const std::exception &e) {
    fmt::print(stderr, "Error: {}\n", e.what());
  }
  return 1;
}



extern "C" int main2(int, char** argv, CCallbacks cb = {}) {
  // Solver should be destroyed after any IloException is handled.
  try {
    typedef mp::SolverApp<mp::IlogCPSolver> IlogCPApp;
    std::unique_ptr<IlogCPApp> s;
    try {
      s.reset(new IlogCPApp());
      s->GetCallbacks() = cb;
      return s->Run(argv);
    }
    catch (const IloException& e) {
      fmt::print(stderr, "Error: {}\n", e.getMessage());
    }
  }
  catch (const std::exception& e) {
    fmt::print(stderr, "Error: {}\n", e.what());
  }
  return 1;
}
