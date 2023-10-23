#ifndef NLWRITER2MISC_H
#define NLWRITER2MISC_H

#include <string>
#include <algorithm>
#include <cstdio>
#include <cstdarg>
#include <cassert>

#include "mp/nl-header.h"
#include "mp/nl-utils2.h"

namespace mp {

/// Text formatter
class TextFormatter {
public:
  /// Construct
  TextFormatter(NLUtils& u,
                bool nlc,
                int outprec=0)
    : utils_(u),
      nl_comments(nlc),
      output_prec(outprec)
  { }
  /// Mode query
  NLHeader::Format Mode() const { return NLHeader::TEXT; }
  /// Text printf
  int apr(File&, const char*, ...);
  /// Text nput
  void nput(File&, double);
  /// Retrieve utils
  NLUtils& Utils() { return utils_; }


private:
  NLUtils& utils_;
  const bool nl_comments{false};
  const int output_prec{0};  // 0 means full precision
};


/// Binary formatter
class BinaryFormatter {
public:
  /// Construct
  BinaryFormatter(NLUtils& u, bool, int)
    : utils_(u) { }
  /// Mode query
  NLHeader::Format Mode() const { return NLHeader::BINARY; }
  /// Binary printf
  int apr(File&, const char*, ...);
  /// Binary nput
  void nput(File&, double);
  /// Retrieve utils
  NLUtils& Utils() { return utils_; }


private:
  NLUtils& utils_;
};


class StringFileWriter;


/// NLWriter2 parameters.
/// @param Formatter: low-level writer: binary or text.
/// @param Feeder2: a class implementing the NLFeeder2
/// concept that feeds model information on request.
/// @param Utils: writer utilities
template <class Formatter, class Feeder2>
struct NLWriter2Params {
  using FormatterType=Formatter;
  using FeederType=Feeder2;
};

}  // namespace mp

#endif // NLWRITER2MISC_H
