/*
 Utilities for writing AMPL solver drivers.

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

#include "solvers/util/driver.h"

#include <cstring>

#include "solvers/getstub.h"

namespace ampl {

Problem::Problem() : asl_(reinterpret_cast<ASL_fg*>(ASL_alloc(ASL_read_fg))) {}

Problem::~Problem() {
  ASL_free(reinterpret_cast<ASL**>(&asl_));
}

bool Problem::Read(char **&argv, Option_Info *oi) {
  ASL *asl = reinterpret_cast<ASL*>(asl_);
  char *stub = getstub_ASL(asl, &argv, oi);
  if (!stub) {
    usage_noexit_ASL(oi, 1);
    return false;
  }
  FILE *nl = jac0dim_ASL(asl, stub, static_cast<ftnlen>(std::strlen(stub)));
  asl_->i.Uvx_ = static_cast<real*>(Malloc(num_vars() * sizeof(real)));
  asl_->i.Urhsx_ = static_cast<real*>(Malloc(num_cons() * sizeof(real)));
  efunc *r_ops_int[N_OPS];
  for (int i = 0; i < N_OPS; ++i)
    r_ops_int[i] = reinterpret_cast<efunc*>(i);
  asl_->I.r_ops_ = r_ops_int;
  asl_->p.want_derivs_ = 0;
  fg_read_ASL(asl, nl, ASL_allow_CLP);
  asl_->I.r_ops_ = 0;
  return true;
}

int Driver::GetOptions(char **argv, Option_Info *oi) {
  return getopts_ASL(reinterpret_cast<ASL*>(problem_.asl_), argv, oi);
}
}
