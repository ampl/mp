#include <stdlib.h>

#include "nlsol_ex_c_nl.h"

NLHeader_C CAPI_ex_Header(void* pex_void) {
  CAPIExample* pex = (CAPIExample*)pex_void;
  NLHeader_C hdr;

  hdr.pi.num_vars = pex->n_var;
  hdr.pi.num_algebraic_cons = pex->n_con;
  hdr.pi.num_objs = pex->n_obj;
  hdr.pi.num_ranges = 0;
  hdr.pi.num_eqns = 0;
  hdr.pi.num_logical_cons = 0;

  // Setting some other common data (more may be needed)

  /** Total number of nonlinear constraints. */
  hdr.pi.num_nl_cons = 0;
  hdr.pi.num_nl_objs = 0;
  hdr.pi.num_compl_conds = 0;
  hdr.pi.num_nl_compl_conds = 0;
  hdr.pi.num_compl_dbl_ineqs = 0;
  hdr.pi.num_compl_vars_with_nz_lb = 0;

  /** Number of nonlinear network constraints. */
  hdr.pi.num_nl_net_cons = 0;
  hdr.pi.num_linear_net_cons = 0;

  /**
    Number of nonlinear variables in constraints including nonlinear
    variables in both constraints and objectives.
   */
  hdr.pi.num_nl_vars_in_cons = 0;

  /**
    Number of nonlinear variables in objectives including nonlinear
    variables in both constraints and objectives.
   */
  hdr.pi.num_nl_vars_in_objs = 0;

  /** Number of nonlinear variables in both constraints and objectives. */
  hdr.pi.num_nl_vars_in_both = 0;

  // Miscellaneous
  // -------------

  /** Number of linear network variables (arcs). */
  hdr.pi.num_linear_net_vars = 0;

  /** Number of functions. */
  hdr.pi.num_funcs = 0;

  // Information about discrete variables
  // ------------------------------------

  /** Number of linear binary variables. */
  hdr.pi.num_linear_binary_vars = 0;

  /** Number of linear non-binary integer variables. */
  hdr.pi.num_linear_integer_vars = pex->n_var_int;

  /**
    Number of integer nonlinear variables in both constraints and objectives.
   */
  hdr.pi.num_nl_integer_vars_in_both = 0;

  /** Number of integer nonlinear variables just in constraints. */
  hdr.pi.num_nl_integer_vars_in_cons = 0;

  /** Number of integer nonlinear variables just in objectives. */
  hdr.pi.num_nl_integer_vars_in_objs = 0;

  // Information about nonzeros
  // --------------------------

  /** Number of nonzeros in constraints' Jacobian. */
  hdr.pi.num_con_nonzeros = pex->n_con_nz;

  /** Number of nonzeros in all objective gradients. */
  hdr.pi.num_obj_nonzeros = pex->n_obj_nz;

  // Information about names
  // -----------------------

  /** Length of longest con/obj name if names are available. */
  hdr.pi.max_con_name_len = 0;    // no need to set

  /** Length of longest variable name if names are available. */
  hdr.pi.max_var_name_len = 0;    // no need to set

  // Information about common expressions
  // ------------------------------------

  /**
    Number of common expressions that appear both in constraints
    and objectives.
   */
  hdr.pi.num_common_exprs_in_both = 0;

  /**
    Number of common expressions that appear in multiple constraints
    and don't appear in objectives.
   */
  hdr.pi.num_common_exprs_in_cons = 0;

  /**
    Number of common expressions that appear in multiple objectives
    and don't appear in constraints.
   */
  hdr.pi.num_common_exprs_in_objs = 0;

  /**
    Number of common expressions that only appear in a single constraint
    and don't appear in objectives.
   */
  hdr.pi.num_common_exprs_in_single_cons = 0;

  /**
    Number of common expressions that only appear in a single objective
    and don't appear in constraints.
   */
  hdr.pi.num_common_exprs_in_single_objs = 0;

  hdr.nli.prob_name = "c_api_example_model";

  return hdr;
}

NLFeeder2_C MakeNLFeeder2_C(
    CAPIExample* pex, int binary) {
  NLFeeder2_C result;

  result.p_user_data_ = pex;
  pex->binary_nl = binary;

  result.Header = CAPI_ex_Header;

  return result;
}

void DestroyNLFeeder2_C(NLFeeder2_C* pf) {
  pf->p_user_data_ = NULL;
}
