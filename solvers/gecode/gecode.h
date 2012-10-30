#ifndef AMPL_SOLVERS_GECODE_H
#define AMPL_SOLVERS_GECODE_H

#include <memory>

#include "solvers/util/driver.h"

struct Option_Info;
struct ASL_fg;

namespace ampl {

// The Gecode driver for AMPL.
class GecodeDriver : public Driver {
 private:
  std::auto_ptr<Option_Info> oinfo_;

 public:
  GecodeDriver();

  // Run the driver.
  int run(char **argv);
};

}

#endif // AMPL_SOLVERS_GECODE_H
