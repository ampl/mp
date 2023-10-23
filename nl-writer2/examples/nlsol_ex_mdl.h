#ifndef NLSOL_EX_MDL_H
#define NLSOL_EX_MDL_H

#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <cassert>

#include "mp/nl-header.h"

/**
 * @brief The ExampleModel struct.
 *
 * A non-linear model with defined variables.
 * Illustrates variable/constraint order,
 * handling of nonlinearities,
 * order of defined variables and their splitting
 * into linear and non-linear parts.
 *

## eXample for NLWriter2/SOLReader2.

## Variable bounds.
## Defined variabes.
## Pure linear constraint.
## Quadratic constraint.
## General nonlinear constraint.
## Initial guess (useful for local solvers.)
## Suffix.
##
## To write NL file and name files in AMPL, use commands
##    ampl: option nl_comments 1;
##    ampl: option auxfiles rc;
##    ampl: write gmodel;
##
## Full documentation on NL format:
## technical report "Writing .nl Files"
## (https://ampl.github.io/nlwrite.pdf.)

var x >=1.4, <= 30;
var y >=0;

var t1 = x + 2*y + y^2;
var t2 = 3*x - 5*y;                              ## linear dvar
var t3 = 4*y + abs(y-2) + 6.38;

minimize Obj_t2t3: 5*x + 3*t2 + 8*t3 + 15;

s.t. C0_lin: 300 >= 5*x - 17*y >= -4;            ## range constraint
s.t. C1_t2: x + 12*y - 18*t2 + 5*t2^2 <= 119;    ## quadratic constraint
s.t. C2_t1t2t3: -38.2*y^2 +109*sin(t2) * (t1 + t3^1.5) >= -30;

let x := 1.5;
let y := 0.11;

suffix zork;
let C1_t2.zork := 5;
 *
 */
struct ExampleModel {
	/// Return NL header.
	mp::NLHeader Header() const {
		mp::NLHeader hdr;

		hdr.num_vars = n_var;
		hdr.num_algebraic_cons = n_con;
		hdr.num_objs = n_obj;
		hdr.num_ranges = 1;
		hdr.num_eqns = 0;
		hdr.num_logical_cons = 0;

		/** Total number of nonlinear constraints. */
		hdr.num_nl_cons = 2;
		hdr.num_nl_objs = 1;
		hdr.num_compl_conds = 0;
		hdr.num_nl_compl_conds = 0;
		hdr.num_compl_dbl_ineqs = 0;
		hdr.num_compl_vars_with_nz_lb = 0;

		/** Number of nonlinear network constraints. */
		hdr.num_nl_net_cons = 0;
		hdr.num_linear_net_cons = 0;

		/**
			Number of nonlinear variables in constraints including nonlinear
			variables in both constraints and objectives.
		 */
		hdr.num_nl_vars_in_cons = 2;

		/**
			Number of nonlinear variables in objectives including nonlinear
			variables in both constraints and objectives.
		 */
		hdr.num_nl_vars_in_objs = 1;

		/** Number of nonlinear variables in both constraints and objectives. */
		hdr.num_nl_vars_in_both = 1;

		// Miscellaneous
		// -------------

		/** Number of linear network variables (arcs). */
		hdr.num_linear_net_vars = 0;

		/** Number of functions. */
		hdr.num_funcs = 0;

		// Information about discrete variables
		// ------------------------------------

		/** Number of linear binary variables. */
		hdr.num_linear_binary_vars = 0;

		/** Number of linear non-binary integer variables. */
		hdr.num_linear_integer_vars = 0;

		/**
			Number of integer nonlinear variables in both constraints and objectives.
		 */
		hdr.num_nl_integer_vars_in_both = 0;

		/** Number of integer nonlinear variables just in constraints. */
		hdr.num_nl_integer_vars_in_cons = 0;

		/** Number of integer nonlinear variables just in objectives. */
		hdr.num_nl_integer_vars_in_objs = 0;

		// Information about nonzeros
		// --------------------------

		/** Number of nonzeros in constraints' Jacobian. */
		hdr.num_con_nonzeros = n_con_nz;

		/** Number of nonzeros in all objective gradients. */
		hdr.num_obj_nonzeros = n_obj_nz;

		// Information about names
		// -----------------------

		/** Length of longest con/obj name if names are available. */
		hdr.max_con_name_len = 0;    // no need to set

		/** Length of longest variable name if names are available. */
		hdr.max_var_name_len = 0;    // no need to set

		// Information about common expressions
		// ------------------------------------

		/**
			Number of common expressions that appear both in constraints
			and objectives.
		 */
		hdr.num_common_exprs_in_both = 1;

		/**
			Number of common expressions that appear in multiple constraints
			and don't appear in objectives.
		 */
		hdr.num_common_exprs_in_cons = 1;

		/**
			Number of common expressions that appear in multiple objectives
			and don't appear in constraints.
		 */
		hdr.num_common_exprs_in_objs = 0;

		/**
			Number of common expressions that only appear in a single constraint
			and don't appear in objectives.
		 */
		hdr.num_common_exprs_in_single_cons = 2;

		/**
			Number of common expressions that only appear in a single objective
			and don't appear in constraints.
		 */
		hdr.num_common_exprs_in_single_objs = 0;

		hdr.prob_name = "example_model";

		return hdr;
	}

  static constexpr int n_var = 2;
  static constexpr int n_con = 3;
  static constexpr int n_obj = 1;

  /////////////////////////////////////////////////
  //////////// Bounds and linear parts ////////////
  /////////////////////////////////////////////////

  /// Variables.
  /// Put y first because it participates non-linearly
  /// both in constraints and the objective.
  /// This might be assumed by some NL solvers.
  double var_lb[n_var] = {       0, 1.4};
  double var_ub[n_var] = {INFINITY,  30};
  const char* var_name[n_var] = {"y", "x"};

  /// Defined variables.
  /// It proves essential (IPOPT, MINOS) to split
  /// those defined variables, whose linear parts
  /// can be substituted into the constraints.
  /// In this example, this is done for t3: the non-linear part
  /// is extracted into a new defined variable nl(t3).
  /// Note that x, y are swapped!
  static constexpr int n_dvar = 4;
  double dvar_linpart[n_dvar][n_var] = {
    {0, 0}, {-5, 3}, {4, 0}, {2, 1}
  };
  const char* dvar_name[n_dvar]
    = {"nl(t3)", "t2", "t3", "t1"};

  /// For some non-linear solvers (MINOS) it is necessary
  /// to provide defined variables, which are used
  /// in a single obj/con, just before that.
  /// This info is returned here:
  /// @param i: i<0 for obj -i-1,
  ///   i>0 for con i-1,
  ///   i=0 for all other defvars.
  /// @return list of defined variables for \a i.
  std::vector<int> DefVarsInItem(int i) const {
    switch (i) {
    case -1: return {};
    case 0: return {0, 1};
    case 1:
    case 3: return {};
    case 2: return {2, 3};
    default:
      assert(0);
      return {};
    }
  }

  /// Algebraic constraint bounds.
  /// Constraints are reordered: nonlinear first.
  /// If we had logical constraints, they'd go last
  /// (represented by expressions only.)
  double con_lb[n_con] = {-INFINITY,      -30,  -4};
  double con_ub[n_con] = {      119, INFINITY, 300};
  const char* con_name[n_con]
                       = {"C1_t2", "C2_t1t2t3", "C0_lin"};

  /// Typedef SparseVec<>
  template <class Elem>
  using SparseVec = std::vector< std::pair<int, Elem> >;

  /// Algebraic constraints: linear parts.
  ///
  /// Linear parts (columns of the Jacobian)
  /// are presented sparse, but should include entries
  /// which can become nonzero in the Jacobian.
  /// Thus, for C2_t1t2t3 we transmit the 0 elements.
  ///
  /// In the linear part of C1_t2 we substitute t2.
  /// We could have put it into the nonlinear expression.
  SparseVec<double> con_linpart[n_con] = {
    { {0, 102}, {1, -53} },
    { {0, 0}, {1, 0} },
    { {0, -17}, {1, 5} }
  };
  /// Sizes of all Jacobian columns except the last
  int col_sizes[n_var-1] = {3};
  int n_con_nz = 6;

  /// Obj sense: min
  int obj_sense = 0;
  /// Objective: linear part.
  /// Similar to the Jacobian,
  /// need all elements which can become nonzero.
  ///
  /// We substitute t2 and the linear part of t3 here.
  SparseVec<double> obj_linpart
      = { {0, 17}, {1, 14} };
  const char* obj_name = "Obj_t2t3";
  int n_obj_nz = 2;

  /////////////////////////////////////////////////
  ////// Non-linear parts: expression trees ///////
  /////////////////////////////////////////////////

  // Here we write expression trees "on the fly".
  // In the case when we'd write expressions from
  // an intermediately stored tree,
  // we might use ExprWriter::EPut(expr_node)
  // to process subtrees recursively.

  /// Write dvar expressions.
  template <class EWriter>
  void WriteDVarExpr(int i, EWriter& ew) const {
    switch (i) {
    case 0:                  // nl(t3)
    {
      auto ew01 = ew.OPut1(15, "abs");
      {
        auto ew0101 = ew01.OPut2(0, "+");
        ew0101.NPut(-2);
        ew0101.VPut(0, "y");
      }
    } break;
    case 1:                  // t2
      ew.NPut(0);            // no non-linear part
      break;
    case 2:                  // t3
    {
      auto ew01 = ew.OPut2(0, "+");
      ew01.VPut(2, "nl(t3)");
      ew01.NPut(6.38);
    } break;
    case 3:                  // t1
    {
      auto ew01 = ew.OPut2(5, "^");   // binary opcode
      ew01.VPut(0, "y");     // 1st arg
      ew01.NPut(2);          // 2nd arg
    } break;
    default:
      assert(false);
    {}
    }
  }

  /// Write constraint non-linear expressions.
  template <class EWriter>
  void WriteConExpr(int i, EWriter& ew) const {
    switch (i) {
    case 0:                  // C1_t2
    {
      auto ew01 = ew.OPut2(2, "*");
      ew01.NPut(5);
      {
        auto ew0101 = ew01.OPut2(5, "^");
        ew0101.VPut(3, "t2");
        ew0101.NPut(2);
      }
    } break;
    case 1:                  // C2_t1t2t3
    {
      auto ew01 = ew.OPut2(0, "+");
      {
        auto ew0101 = ew01.OPut2(2, "*");
        ew0101.NPut(-38.2);
        {
          auto ew010101 = ew0101.OPut2(5, "^");
          ew010101.VPut(0, "y");
          ew010101.NPut(2);
        }
      }
      {
        auto ew0102 = ew01.OPut2(2, "*");
        {
          auto ew010201 = ew0102.OPut2(2, "*");
          ew010201.NPut(109);
          {
            auto ew01020101 = ew010201.OPut1(41, "sin");
            ew01020101.VPut(3, "t2");
          }
        }
        {
          auto ew010202 = ew0102.OPut2(0, "+");
          ew010202.VPut(5, "t1");
          {
            auto ew01020201 = ew010202.OPut2(5, "^");
            ew01020201.VPut(4, "t3");
            ew01020201.NPut(1.5);
          }
        }
      }
    } break;
    case 2:                  // C0_lin
      ew.NPut(0);            // no non-linear part
      break;
    default:
      assert(false);
    {}
    }
  }

  template <class EWriter>
  void WriteObjExpr(EWriter& ew) const {
    auto ew01 = ew.OPut2(0, "+");
    {
      auto ew0101 = ew01.OPut2(2, "*");
      ew0101.NPut(8);
      ew0101.VPut(2, "nl(t3)");
    }
    ew01.NPut(66.04);
  }

  /// Primal initial guess.
  double ini_x[n_var] = {0.11, 1.5};

  /// Suffixes.
  struct Suffix {
    const char* name_;
    int kind_;
    std::vector<double> values_;   // store double's always
  };
  static constexpr int n_suf = 1;
  Suffix suf[n_suf] = { {
                          "zork",
                          1,          // constraints, integer
                          {5, 0, 0}
                        } };


  /////////////////////////////////////////////////
  /////////////////// Solution ////////////////////
  /////////////////////////////////////////////////

  std::vector<double> sol_y_;
  std::vector<double> sol_x_;
  int objno_ = -2;
  int solve_result_ = -1;
  std::map<
    std::pair< std::string, int >,
    std::vector< double >
  > suf_out_;
};

#endif // NLSOL_EX_MDL_H
