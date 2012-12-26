/*
 Implementation of the Args class.

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

#include "tests/args.h"

#include <cstring>

void Args::Add(const char *arg) {
  if (!arg) return;
  ++argc_;
  store_.insert(store_.end(), arg, arg + std::strlen(arg) + 1);
}

Args::Args(const char *arg1, const char *arg2,
      const char *arg3, const char *arg4)
: argc_(0) {
  Add(arg1);
  Add(arg2);
  Add(arg3);
  Add(arg4);
}

Args::operator char **() {
  argv_.resize(argc_ + 1);
  for (std::size_t i = 0, j = 0; i < argc_;
      j += std::strlen(&store_[j]) + 1, ++i) {
    argv_[i] = &store_[j];
  }
  return &argv_[0];
}
