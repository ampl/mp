#ifndef PROBLEM_FLATTENER_BASE_H
#define PROBLEM_FLATTENER_BASE_H

namespace mp {

/// An abstract base for ProblemFlattener
class BasicProblemFlattener {
public:
  /// Destructor
  virtual ~BasicProblemFlattener() { }

  /// Number of variables in the original model
  virtual int num_vars_orig() const = 0;
};

}  // namespace mp

#endif // PROBLEM_FLATTENER_BASE_H
