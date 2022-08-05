#ifndef BACKENDWITHVALCVT_H
#define BACKENDWITHVALCVT_H

/**
 * Backends using ValuePresolver should derive from here
 */

#include "mp/valcvt-base.h"

namespace mp {

class BackendWithValuePresolver :
    public BasicValuePresolverKeeper {
public:
  // otherwise empty
};

}  // namespace mp

#endif // BACKENDWITHVALCVT_H
