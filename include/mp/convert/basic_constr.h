#ifndef BASIC_CONSTR_H
#define BASIC_CONSTR_H

#include <vector>


namespace mp {

/// Custom constraints to derive from
class BasicConstraint {
public:
};

/// A constraint whose arguments are an array of variables
class VarArrayArgConstraint : public BasicConstraint {
public:
  using ArgArray = std::vector<int>;       // by default, our arguments are an array of variables
private:
  ArgArray args_;
public:
  VarArrayArgConstraint(ArgArray&& aa) : args_(aa) {}
  const ArgArray& GetArguments() const { return args_; }
};

/// A constraint extension which defines a variable
class DefiningConstraint {
  int result_var_=-1;                // defined var can be optional
public:
  DefiningConstraint(int v) : result_var_(v) {}
  int GetResultVar() const { return result_var_; }
};

} // namespace mp

#endif // BASIC_CONSTR_H
