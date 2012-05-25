#ifndef AMPL_SOLVERS_CONCERT_H
#define AMPL_SOLVERS_CONCERT_H

#include <stdexcept>

struct expr;

class IloConstraint;
class IloExpr;
class IloEnv;
class IloModel;
class IloNumVar;
class IloNumVarArray;

extern IloEnv env;
extern IloModel mod;
extern IloNumVarArray Var;

extern int usenumberof;
extern int debugexpr;

// Variable subscripted by a variable - not implemented in AMPL yet.
enum { OPVARSUBVAR = 99 };

// Operation types
enum {
  // Unary operation
  OPTYPE_UNARY = 1,

  // Binary operation
  OPTYPE_BINARY = 2,

  // Variable-argument function such as min or max
  OPTYPE_VARARG = 3,

  // Piecewise-linear term
  OPTYPE_PLTERM = 4,

  // The if-then-else expression
  OPTYPE_IF = 5,

  // The sum expression
  OPTYPE_SUM = 6,

  // Function call
  OPTYPE_FUNCALL = 7,

  // String
  OPTYPE_STRING = 8,

  // Number
  OPTYPE_NUMBER = 9,

  // Variable
  OPTYPE_VARIABLE = 10,

  // The count expression
  OPTYPE_COUNT = 11
};

IloExpr build_expr(expr *e);
IloConstraint build_constr (expr*);
IloNumVar build_numberof (expr*);
bool same_expr (expr *e1, expr *e2);

const char *get_opname(int opcode);

int concert_main(int argc, char **argv);

// Thrown when an ASL expression cannot be converted to a Concert one.
class UnsupportedExprError : public std::runtime_error {
 public:
  UnsupportedExprError(const char *expr) :
    std::runtime_error(std::string("unsupported expression: ") + expr) {}
};


#endif // AMPL_SOLVERS_CONCERT_H
