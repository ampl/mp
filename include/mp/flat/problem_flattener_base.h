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

  /// Original variable's lower bound
  virtual double var_orig_lb(int i) const = 0;
  /// Original variable's upper bound
  virtual double var_orig_ub(int i) const = 0;

  /// Quadratize ^2?
  virtual bool IfQuadratizePow2() const = 0;

  /// Mutliply-out cardinality
  virtual double MultOutCard() const = 0;
};

}  // namespace mp

#endif // PROBLEM_FLATTENER_BASE_H
