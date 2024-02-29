/*
 JSON write utils.

 Copyright (C) 2024 AMPL Optimization Inc.

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

 Author: Gleb Belov
 */

#ifndef UTILJSONWRITE_HPP
#define UTILJSONWRITE_HPP

#include <cassert>

#include "mp/util-json-write.h"

namespace mp {

template <class StringWriter>
typename MiniJSONWriter<StringWriter>::Node
MiniJSONWriter<StringWriter>::operator++() {
  EnsureArray();
  InsertElementSeparator();
  return Node{*this};
}

template <class StringWriter>
typename MiniJSONWriter<StringWriter>::Node
MiniJSONWriter<StringWriter>::operator[](const char* key) {
  EnsureDictionary();
  InsertElementSeparator();
  wrt_.write("{}: ",  key);
  return Node{*this};
}

template <class StringWriter>
void MiniJSONWriter<StringWriter>::EnsureUnset() {
  assert(Kind::Unset == kind_);
}

template <class StringWriter>
void MiniJSONWriter<StringWriter>::EnsureScalar() {
  EnsureUnset();
  kind_ = Kind::Scalar;
}

template <class StringWriter>
void MiniJSONWriter<StringWriter>::EnsureArray() {
  if (Kind::Unset == kind_) {
    kind_ = Kind::Array;
    wrt_.write("[");
  } else {
    assert(Kind::Array == kind_);    // cannot change
  }
}

template <class StringWriter>
void MiniJSONWriter<StringWriter>::EnsureDictionary() {
  if (Kind::Unset == kind_) {
    kind_ = Kind::Dict;
    wrt_.write("{}", '{');
  } else {
    assert(Kind::Dict == kind_);    // cannot change
  }
}

template <class StringWriter>
void MiniJSONWriter<StringWriter>::EnsureCanWrite() {
  switch (kind_) {
  case Kind::Unset:
    assert(false && "unset MiniJSONWriter node kind for writing");
    break;
  case Kind::Scalar:
    assert(0==n_written_);
    break;
  case Kind::Array:
  case Kind::Dict:
    break;
  case Kind::Closed:     // done already
    assert(false && "writing into closed MiniJSONWriter node");
    break;
  default:
    assert(false && "unknown MiniJSONWriter node kind");
    break;
  }
}

template <class StringWriter>
void MiniJSONWriter<StringWriter>::InsertElementSeparator() {
  if (n_written_)
    wrt_.write(", ");
}

template <class StringWriter>
void MiniJSONWriter<StringWriter>::Close() {
  switch (kind_) {
  case Kind::Unset:
    wrt_.write("[]");  // empty array. Is this ok to default to?
    break;
  case Kind::Scalar:
    assert(1==n_written_);
    break;
  case Kind::Array:
    wrt_.write("]");
    break;
  case Kind::Dict:
    wrt_.write("{}", '}');
    break;
  case Kind::Closed:     // done already
    break;
  default:
    assert(false && "unknown MiniJSONWriter node kind");
    break;
  }
  kind_ = Kind::Closed;
}


}  // namespace mp

#endif // UTILJSONWRITE_HPP
