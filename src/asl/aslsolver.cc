/*
 A mathematical optimization solver.

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

extern "C" {
#include "solvers/getstub.h"
}

#undef Char
#undef filename
#undef ampl_vbtol
#undef write_sol

#include "mp/clock.h"
#include "aslproblem.h"
#include "aslsolver.h"

mp::ASLSolver::ASLSolver(
    fmt::StringRef name, fmt::StringRef long_name, long date, int flags)
  : SolverImpl<internal::ASLBuilder>(name, long_name, date, flags) {
}

void mp::ASLSolver::RegisterSuffixes(ASL *asl) {
  std::size_t num_suffixes = suffixes_.size();
  std::vector<SufDecl> suffix_decls(num_suffixes);
  for (std::size_t i = 0; i < num_suffixes; ++i) {
    const SuffixInfo &si = suffixes_[i];
    SufDecl &sd = suffix_decls[i];
    sd.name = const_cast<char*>(si.name());
    sd.table = const_cast<char*>(si.table());
    sd.kind = si.kind();
    sd.nextra = si.nextra();
  }
  if (asl->i.nsuffixes == 0 && num_suffixes != 0)
    suf_declare_ASL(asl, &suffix_decls[0], static_cast<int>(num_suffixes));
}

mp::Problem::Proxy mp::ASLSolver::GetProblemBuilder(fmt::StringRef stub) {
  Problem::Proxy proxy(ASL_alloc(ASL_read_fg), read_flags_);
  std::size_t stub_len = stub.size();
  Edaginfo &info = proxy.asl_->i;
  info.filename_ = reinterpret_cast<char*>(M1alloc_ASL(&info, stub_len + 5));
  std::strcpy(info.filename_, stub.c_str());
  info.stub_end_ = info.filename_ + stub_len;
  std::strcpy(info.filename_ + stub_len, ".nl");
  RegisterSuffixes(proxy.asl_);
  return proxy;
}

void mp::ASLSolver::Solve(Problem &p, SolutionHandler &sh) {
  RegisterSuffixes(p.asl_);
  ASLSolutionHandler asl_sol_handler(sh, p);
  DoSolve(p, asl_sol_handler);
}
