/*
 AMPL suffix support
 Suffixes are values associated with model components.
 See http://www.ampl.com/NEW/suffixes.html

 Copyright (C) 2015 AMPL Optimization Inc

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use	, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#include "mp/suffix.h"

mp::SuffixSet::~SuffixSet() {
  // Deallocate names and values.
  for (Set::iterator i = set_.begin(), e = set_.end(); i != e; ++i) {
    delete [] i->name.c_str();
    if ((i->kind & suf::FLOAT) != 0)
      delete [] i->dbl_values;
    else
      delete [] i->int_values;
  }
}

mp::Suffix::Impl *mp::SuffixSet::DoAdd(fmt::StringRef name,
                                       int kind, int num_values) {
  Suffix::Impl *impl = const_cast<Suffix::Impl*>(
        &*set_.insert(Suffix::Impl(name, kind)).first);
  // Set name to empty string so that it is not deleted if new throws.
  std::size_t size = name.size();
  impl->name = fmt::StringRef(0, 0);
  char *name_copy = new char[size + 1];
  const char *s = name.c_str();
  std::copy(s, s + size, fmt::internal::make_ptr(name_copy, size));
  name_copy[size] = 0;
  impl->name = fmt::StringRef(name_copy, size);
  impl->num_values = num_values;
  return impl;
}
