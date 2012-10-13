#ifndef AMPL_SOLVERS_GECODE_H
#define AMPL_SOLVERS_GECODE_H

#include <memory>

struct Option_Info;
struct ASL_fg;

// The Gecode driver for AMPL.
class Driver {
 private:
  ASL_fg *asl;
  std::auto_ptr<Option_Info> oinfo_;

  // Do not implement.
  Driver(const Driver&);
  Driver& operator=(const Driver&);

 public:
  Driver();
  virtual ~Driver();

  // Run the driver.
  int run(char **argv);
};

#endif // AMPL_SOLVERS_GECODE_H
