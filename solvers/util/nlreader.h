/*
 An nl reader.

 Copyright (C) 2013 AMPL Optimization Inc

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#include "solvers/util/error.h"

namespace ampl {

class ParseError : public Error {
 private:
  std::string filename_;
  int line_;
  int column_;

 public:
  ParseError(fmt::StringRef filename,
      int line, int column, fmt::StringRef message)
  : Error(message), filename_(filename), line_(line), column_(column) {}
  ~ParseError() FMT_NOEXCEPT(true) {}

  const std::string &filename() const { return filename_; }
  int line() const { return line_; }
  int column() const { return column_; }
};

class TextReader;

enum {
  MAX_NL_OPTIONS = 9,
  VBTOL_OPTION   = 1,
  READ_VBTOL     = 3
};

// NL file header.
struct NLHeader {
  int num_options;
  int options[MAX_NL_OPTIONS];

  // Extra info for writing solution.
  double ampl_vbtol;

  // Total number of variables.
  int num_vars;

  // Total number of constraints including ranges, equality constraints
  // and logical constraints.
  int num_cons;

  // Total number of objectives.
  int num_objs;

  // Number of ranges (constraints with -Infinity < LHS < RHS < Infinity).
  int num_ranges;

  // Number of equality constraints or -1 if unknown (AMPL prior to 19970627).
  int num_eqns;

  // Number of logical constraints.
  int num_logical_cons;

  // Nonlinear and complementarity information
  // -----------------------------------------

  // Total number of nonlinear constraints.
  int num_nl_cons;

  // Total number of nonlinear objectives.
  int num_nl_objs;

  // Total number of complementarity conditions.
  int num_compl_conds;

  // Number of nonlinear complementarity conditions.
  int num_nl_compl_conds;

  // Number of complementarities involving double inequalities
  // (for ASL_cc_simplify).
  int num_compl_dbl_ineqs;

  // Number of complemented variables with a nonzero lower bound
  // (for ASL_cc_simplify).
  int num_compl_vars_with_nz_lb;

  // Information about network constraints
  // -------------------------------------

  // Number of nonlinear network constraints.
  int num_nl_net_cons;

  // Number of linear network constraints.
  int num_linear_net_cons;

  // Information about nonlinear variables
  // -------------------------------------

  // Number of nonlinear variables in constraints including nonlinear
  // variables in both constraints and objectives.
  int num_nl_vars_in_cons;

  // Number of nonlinear variables in objectives including nonlinear
  // variables in both constraints and objectives.
  int num_nl_vars_in_objs;

  // Number of nonlinear variables in both  constraints and objectives.
  int num_nl_vars_in_both;

  // Miscellaneous
  // -------------

  // Number of linear network variables (arcs).
  int num_linear_net_vars;

  // Number of functions.
  int num_funcs;

  // Flags: 1 = want output suffixes.
  int flags;

  // Information about discrete variables
  // ------------------------------------

  // Number of linear binary variables.
  int num_linear_binary_vars;

  // Number of linear non-binary integer variables.
  int num_linear_integer_vars;

  // Number of integer nonlinear variables in both constraints and objectives.
  int num_nl_integer_vars_in_both;

  // Number of integer nonlinear variables just in constraints.
  // TODO: make consistent with num_nl_vars_in_cons
  int num_nl_integer_vars_in_cons;

  // Number of integer nonlinear variables just in objectives.
  // TODO: make consistent with num_nl_vars_in_objs
  int num_nl_integer_vars_in_objs;

  // Information about nonzeros
  // --------------------------

  // Number of nonzeros in constraints' Jacobian.
  int num_con_nonzeros;


  // Number of nonzeros in all objective gradients.
  int num_obj_nonzeros;

  // Information about names
  // -----------------------

  // Length of longest constraint name (if stub.row exists).
  int max_con_name_len;

  // Length of longest variable name (if stub.col exists).
  int max_var_name_len;

  // Information about common expressions
  // ------------------------------------

  // TODO: improve naming
  int num_common_b_exprs;
  int num_common_con_exprs;
  int num_common_obj_exprs;
  int num_common_con1_exprs;
  int num_common_obj1_exprs;
};

class NLHandler {
 public:
  virtual ~NLHandler();
  virtual void HandleHeader(const NLHeader &h) = 0;
};

// An nl reader.
class NLReader {
 private:
  NLHandler *handler_;

  // Reads an expression.
  void ReadExpr(TextReader &reader);

  // Reads a linear expression.
  void ReadLinearExpr(TextReader &reader, int num_terms);

  // Reads bounds.
  void ReadBounds(TextReader &reader, int num_bounds);

  // Read the column offsets, the cumulative sums of the numbers of
  // nonzeros in the first num_var − 1 columns of the Jacobian matrix.
  void ReadColumnOffsets(TextReader &reader, int num_vars);

 public:
  explicit NLReader(NLHandler *h = 0) : handler_(h) {}

  // Reads a problem from a file.
  void ReadFile(fmt::StringRef filename);

  // Reads a problem from a string.
  // name: Name to be used when reporting errors.
  void ReadString(fmt::StringRef str, fmt::StringRef name = "(input)");
};
}
