/**
 * "Implementation" of NLFeeder2_C:
 * A C++ class "wrapping" it
 * in order to interface it for NLWriter2
 */
#ifndef NLFEEDER2CIMPL_H
#define NLFEEDER2CIMPL_H

#include "api/c/nl-feeder2-c.h"

#include "mp/nl-feeder2.h"

namespace mp {

/// Wrap NLFeeder2_C into a C++ class,
/// in order to interface it for NLWriter2
class NLFeeder2_C_Impl
    : public NLFeeder2<NLFeeder2_C_Impl, void*> {
public:
  /// Construct
  NLFeeder2_C_Impl(NLFeeder2_C* pnlf2)
    : p_nlf2_c_(pnlf2) { }

private:
  NLFeeder2_C* p_nlf2_c_ {nullptr};
};

}  // namespace mp

#endif // NLFEEDER2CIMPL_H
