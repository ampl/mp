#ifndef BASIC_CONSTR_H
#define BASIC_CONSTR_H

#include <vector>


namespace mp {

/// Custom constraints to derive from, so that overloaded default settings work
class BasicConstraint {
public:
};

/// A constraint whose arguments are an array of variables
class VarArrayArgConstraint : public BasicConstraint {
public:
  using Arguments = std::vector<int>;       // by default, our arguments are an array of variables
private:
  Arguments args_;
public:
  VarArrayArgConstraint() { }
  VarArrayArgConstraint(Arguments&& aa) : args_(std::move(aa)) {}
  const Arguments& GetArguments() const { return args_; }
  Arguments& GetArguments() { return args_; }
};

/// A constraint extension which defines a variable
class DefiningConstraint {
  int result_var_=-1;                // defined var can be optional
public:
  DefiningConstraint(int v=-1) : result_var_(v) {}
  void SetResultVar(int v) { result_var_=v; }
  int GetResultVar() const { return result_var_; }
};

} // namespace mp

#endif // BASIC_CONSTR_H
