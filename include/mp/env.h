#ifndef ENV_H
#define ENV_H

#include "mp/solver-base.h"

namespace mp {

/// Typedef Env to access solver options, output handler, etc
using Env = BasicSolver;

/// Facilitate storing a reference to Env
/// to access solver options, output handler, etc
class EnvKeeper {
public:
  /// Construct
  EnvKeeper(Env& e) : env_(e) { }

  /// GetEnv() const
  const Env& GetEnv() const { return env_; }
  /// GetEnv()
  Env& GetEnv() { return env_; }


private:
  Env& env_;
};


} // namespace mp

#endif // ENV_H
