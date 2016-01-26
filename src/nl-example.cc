// NL reader example

#include "mp/nl.h"
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>

namespace expr = mp::expr;

enum Expr {OTHER, CONST};

// Count the number of divisions by non-constant expressions.
struct ExprCounter : mp::NLHandler<Expr> {
  int num_divs;
  ExprCounter() : num_divs(0) {}
  Expr OnBinary(mp::expr::Kind kind, Expr, Expr rhs) {
    if (kind == mp::expr::DIV && rhs != CONST)
      ++num_divs;
    return OTHER;
  }
  Expr OnNumericConstant(double) { return CONST; }
};

// Print problem dimensions: the number of constraints, variables and nonzeros
// (in constraints).
struct DimensionPrinter : mp::NLHandler<int> {
  void OnHeader(const mp::NLHeader &h) {
    fmt::print("Variables:   {}\n", h.num_vars);
    fmt::print("Constraints: {}\n", h.num_algebraic_cons);
    fmt::print("Nonzeros:    {}\n", h.num_con_nonzeros);
  }
};

// Print sparsity pattern of the objective and constraint gradients.
struct SparsityPrinter : mp::NLHandler<int> {
  struct LinearObjHandler {
    void AddTerm(int var_index, double) {
      fmt::print("  {}\n", var_index);
    }
  };
  typedef LinearObjHandler LinearConHandler;

  LinearObjHandler OnLinearObjExpr(int obj_index, int) {
    fmt::print("Objective {}:\n", obj_index);
    return LinearObjHandler();
  }

  LinearConHandler OnLinearConExpr(int con_index, int) {
    fmt::print("Constraint {}:\n", con_index);
    return LinearConHandler();
  }
};

// Print objective or constraint expression in SSA-like form.
class ExprPrinter : public mp::NLHandler<std::string> {
 private:
  // Next expression ID for SSA-like output.
  int next_expr_id_;

  // Output for the currectly formatted expression.
  fmt::MemoryWriter output_;

  // The formatted nonlinear expression ID.
  int nonlinear_expr_id_;

  // Index of the objective to format or -1.
  int obj_index_;

  // Index of the constraint to format or -1.
  int con_index_;

  // Formats arguments and appends the result to ``output``.
  template <typename... T>
  std::string format(const char *format_str, const T &... args) {
    std::string id = fmt::format("e{}", next_expr_id_++);
    output_ << id << " = ";
    output_.write(format_str, args...);
    output_ << '\n';
    return id;
  }

  // Reset the writer for the next expression.
  void reset() {
    output_.clear();
    next_expr_id_ = 0;
  }

 public:
  ExprPrinter(int obj_index, int con_index)
    : next_expr_id_(0), obj_index_(obj_index), con_index_(con_index) {}

  std::string OnNumber(double value) {
    return format("{}", value);
  }

  std::string OnVariableRef(int var_index) {
    return fmt::format("x{}", var_index);
  }

  std::string OnUnary(expr::Kind kind, std::string arg) {
    return kind == expr::MINUS ?
          format("-{}", arg) : format("{}({})", expr::str(kind), arg);
  }

  std::string OnBinary(expr::Kind kind, std::string lhs, std::string rhs) {
    const char *op = expr::str(kind);
    return kind < expr::ATAN2 ?
          format("{} {} {}", lhs, op, rhs) : format("{}({}, {})", op, lhs, rhs);
  }

  struct NumericArgHandler {
    std::string args;
    void AddArg(std::string arg) {
      args = args.empty() ? arg : args + ", " + arg;
    }
  };

  NumericArgHandler BeginSum(int) {
    return NumericArgHandler();
  }

  std::string EndSum(NumericArgHandler handler) {
    return format("sum({})", handler.args);
  }

  void OnObj(int obj_index, mp::obj::Type, std::string) {
    if (obj_index == obj_index_) {
      fmt::print("{}", output_.c_str());
      nonlinear_expr_id_ = next_expr_id_ - 1;
    }
    reset();
  }

  void OnAlgebraicCon(int con_index, std::string) {
    if (con_index == con_index_) {
      fmt::print("{}", output_.c_str());
      nonlinear_expr_id_ = next_expr_id_ - 1;
    }
    reset();
  }

  class LinearPartHandler {
   private:
    int num_terms;  // The remaining number of terms to format.
    int nonlinear_expr_id;
    int obj_index;
    char type; // Type: 'o' - objective, 'c' - constraint.

   public:
    explicit LinearPartHandler(int num_terms = 0, int nonlinear_expr_id = 0,
                               int obj_index = 0, char type = 'o')
      : num_terms(num_terms), nonlinear_expr_id(nonlinear_expr_id),
        obj_index(obj_index), type(type) {}

    void AddTerm(int var_index, double coef) {
      if (num_terms == 0) return;
      --num_terms;
      if (coef != 0) {
        ++nonlinear_expr_id;
        fmt::print("e{} = ", nonlinear_expr_id);
        if (nonlinear_expr_id != 0)
          fmt::print("e{} + ", nonlinear_expr_id - 1);
        fmt::print("x{} * {}\n", var_index, coef);
      }
      if (num_terms == 0)
        fmt::print("{}{} = e{}\n", type, obj_index, nonlinear_expr_id);
    }
  };

  typedef LinearPartHandler LinearObjHandler;

  LinearObjHandler OnLinearObjExpr(int obj_index, int num_linear_terms) {
    if (obj_index != obj_index_)
      return LinearObjHandler();
    if (num_linear_terms == 0)
      fmt::print("o{} = e{}\n", obj_index, nonlinear_expr_id_);
    return LinearObjHandler(num_linear_terms, nonlinear_expr_id_,
                            obj_index, 'o');
  }

  typedef LinearPartHandler LinearConHandler;

  LinearConHandler OnLinearConExpr(int con_index, int num_linear_terms) {
    if (con_index != con_index_)
      return LinearConHandler();
    if (num_linear_terms == 0)
      fmt::print("c{} = e{}\n", con_index, nonlinear_expr_id_);
    return LinearConHandler(num_linear_terms, nonlinear_expr_id_,
                            con_index, 'c');
  }
};

int PrintUsage(const char *executable_name) {
  fmt::print("Usage:\n");
  fmt::print("{} (count|dimen|sparsity) <filename>\n", executable_name);
  fmt::print("  count: count the number of divisions by non-constant "
             "expressions\n");
  fmt::print("  dimen: print problem dimensions\n");
  fmt::print("  sparsity: print sparsity pattern\n");
  fmt::print("{} (obj-expr|con-expr) <filename> <index>\n", executable_name);
  fmt::print("  obj-expr: print objective expression\n");
  fmt::print("  con-expr: print constraint expression\n");
  return EXIT_FAILURE;
}

int main(int argc, char **argv) {
  const char *executable_name = argv[0];
  if (argc < 3)
    return PrintUsage(executable_name);
  const char *command = argv[1], *filename = argv[2];
  if (std::strcmp(command, "count") == 0 && argc == 3) {
    ExprCounter counter;
    mp::ReadNLFile(filename, counter);
    fmt::print("{} has {} division(s) by non-const\n",
               filename, counter.num_divs);
  } else if (std::strcmp(command, "dimen") == 0 && argc == 3) {
    DimensionPrinter printer;
    mp::ReadNLFile(filename, printer);
  } else if (std::strcmp(command, "sparsity") == 0 && argc == 3) {
    SparsityPrinter printer;
    mp::ReadNLFile(filename, printer);
  } else if (std::strcmp(command, "obj-expr") == 0 && argc == 4) {
    ExprPrinter printer(std::atoi(argv[3]), -1);
    mp::ReadNLFile(filename, printer);
  } else if (std::strcmp(command, "con-expr") == 0 && argc == 4) {
    ExprPrinter printer(-1, std::atoi(argv[3]));
    mp::ReadNLFile(filename, printer);
  } else {
    return PrintUsage(executable_name);
  }
}
