#ifndef AMPL_SOLVERS_CONCERT_H
#define AMPL_SOLVERS_CONCERT_H

#include <memory>
#include <vector>
#include <ilconcert/ilomodel.h>

struct expr;
struct keyword;
struct Option_Info;

extern int usenumberof;
extern int debugexpr;

// Variable subscripted by a variable - not implemented in AMPL yet.
enum { OPVARSUBVAR = 99 };

class NumberOf {
 private:
  IloIntVarArray cards_;
  IloIntArray values_;
  IloIntVarArray vars_;
  const expr *numberofexpr_;

 public:
  NumberOf(IloIntVarArray cards, IloIntArray values,
      IloIntVarArray vars, const expr *e) :
    cards_(cards), values_(values), vars_(vars), numberofexpr_(e) {}

  IloInt num_vars() const {
    return vars_.getSize();
  }

  const expr *numberofexpr() const {
    return numberofexpr_;
  }

  IloDistribute to_distribute(IloEnv env) const {
    return IloDistribute (env, cards_, values_, vars_);
  }

  IloIntVar add(double value, IloEnv env);
};

// The Concert driver for AMPL.
class Driver {
 private:
  IloEnv env_;
  IloModel mod_;
  IloNumVarArray vars_;
  IloAlgorithm alg_;
  std::vector<NumberOf> numberofs_;
  std::vector<char> version_;
  std::auto_ptr<Option_Info> oinfo_;
  static keyword keywords_[];

  // Do not implement.
  Driver(const Driver&);
  Driver& operator=(const Driver&);

  // Builds an array of expressions from the argument list of e.
  IloNumExprArray build_minmax_array(const expr *e);

  // Given a node for a number-of operator that has a constant as its first
  // operand, adds it to the driver's data structure that collects these
  // operators.
  IloNumVar build_numberof(const expr *e);

 public:
  // Options accessible from AMPL.
  enum Option {
    DEBUGEXPR,
    ILOGOPTTYPE,
    TIMING,
    USENUMBEROF,
    NUM_OPTIONS
  };

 private:
  int options_[NUM_OPTIONS];

  static char *set_option(Option_Info *oi, keyword *kw, char *value);
  static char *set_option0(Option_Info *oi, keyword *kw, char *value);
  static char *set_option1(Option_Info *oi, keyword *kw, char *value);

 public:
  Driver();
  virtual ~Driver();

  IloEnv env() const { return env_; }
  IloModel mod() const { return mod_; }
  IloAlgorithm alg() const { return alg_; }

  IloNumVarArray vars() const { return vars_; }
  void set_vars(IloNumVarArray vars) { vars_ = vars; }

  int get_option(Option opt) const { return options_[opt]; }
  void use_numberof(bool use = true) { options_[USENUMBEROF] = use; }

  // Walks an expression tree and constructs a function,
  // returning a Concert IloExpr.
  IloExpr build_expr(const expr *e);

  // Walks an expression tree and constructs a constraint,
  // returning a Concert IloConstraint.
  IloConstraint build_constr(const expr *e);

  void finish_building_numberof();

  int run(int argc, char **argv);
};

#endif // AMPL_SOLVERS_CONCERT_H
