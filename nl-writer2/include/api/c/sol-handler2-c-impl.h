/**
 * "Implementation" of SOLHandler2_C:
 * A C++ class "wrapping" it
 * in order to interface it for NLReader2
 */
#ifndef SOLHANDLER2CIMPL_H
#define SOLHANDLER2CIMPL_H

#include "api/c/sol-handler2-c.h"

#include "mp/sol-handler2.h"

namespace mp {

/// Wrap SOLHandler2_C into a C++ class,
/// in order to interface it for NLReader2
class SOLHandler2_C_Impl
    : public SOLHandler2 {
public:
  /// Construct
  SOLHandler2_C_Impl(SOLHandler2_C* psh2)
    : solh2_c_(*psh2) { }

private:
  /// Just store copy
  const SOLHandler2_C solh2_c_;
};

}  // namespace mp

#endif // SOLHANDLER2CIMPL_H
