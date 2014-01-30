/*
 An AMPL expression factory.

 Copyright (C) 2014 AMPL Optimization Inc

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

#include "solvers/util/expr-factory.h"

#include "solvers/util/nlreader.h"

namespace ampl {

ExprFactory::ExprFactory(const NLHeader &h) : asl_(ASL_alloc(ASL_read_fg)) {
  ASL_fg *asl = reinterpret_cast<ASL_fg*>(asl_);
  for (int i = 0; i < N_OPS; ++i)
    r_ops_[i] = reinterpret_cast<efunc*>(i);
  asl->I.r_ops_ = r_ops_;

  asl->i.ampl_options_[0] = h.num_options;
  for (int i = 0; i < ampl::MAX_NL_OPTIONS; ++i)
    asl->i.ampl_options_[i + 1] = h.options[i];
  asl->i.ampl_vbtol_ = h.ampl_vbtol;

  asl->i.n_var_ = h.num_vars;
  asl->i.n_con_ = h.num_algebraic_cons;
  asl->i.n_obj_ = h.num_objs;
  asl->i.nranges_ = h.num_ranges;
  asl->i.n_eqn_ = h.num_eqns;
  asl->i.n_lcon_ = h.num_logical_cons;

  asl->i.nlc_ = h.num_nl_cons;
  asl->i.nlo_ = h.num_nl_objs;
  asl->i.n_cc_ = h.num_compl_conds;
  asl->i.nlcc_ = h.num_nl_compl_conds;
  asl->i.ndcc_ = h.num_compl_dbl_ineqs;
  asl->i.nzlb_ = h.num_compl_vars_with_nz_lb;

  asl->i.nlnc_ = h.num_nl_net_cons;
  asl->i.lnc_ = h.num_linear_net_cons;

  asl->i.nlvc_ = h.num_nl_vars_in_cons;
  asl->i.nlvo_ = h.num_nl_vars_in_objs;
  asl->i.nlvb_ = h.num_nl_vars_in_both;

  asl->i.nwv_ = h.num_linear_net_vars;
  asl->i.nfunc_ = h.num_funcs;
  asl->i.flags = h.flags;

  asl->i.nbv_ = h.num_linear_binary_vars;
  asl->i.niv_ = h.num_linear_integer_vars;
  asl->i.nlvbi_ = h.num_nl_integer_vars_in_both;
  asl->i.nlvci_ = h.num_nl_integer_vars_in_cons;
  asl->i.nlvoi_ = h.num_nl_integer_vars_in_objs;

  asl->i.nzc_ = h.num_con_nonzeros;
  asl->i.nzo_ = h.num_obj_nonzeros;

  asl->i.maxrownamelen_ = h.max_con_name_len;
  asl->i.maxcolnamelen_ = h.max_var_name_len;

  asl->i.comb_ = h.num_common_exprs_in_both;
  asl->i.comc_ = h.num_common_exprs_in_cons;
  asl->i.como_ = h.num_common_exprs_in_objs;
  asl->i.comc1_ = h.num_common_exprs_in_cons1;
  asl->i.como1_ = h.num_common_exprs_in_objs1;

  asl->i.n_var0 = asl->i.n_var1 = asl->i.n_var_;
  asl->i.n_con0 = asl->i.n_con1 = asl->i.n_con_;
  int nlv = asl->i.nlvc_;
  if (nlv < asl->i.nlvo_)
    nlv = asl->i.nlvo_;
  if (nlv <= 0)
    nlv = 1;
  asl->i.x0len_ = nlv * sizeof(double);
  asl->i.x0kind_ = ASL_first_x;
  asl->i.n_conjac_[0] = 0;
  asl->i.n_conjac_[1] = asl->i.n_con_;

  // confusion arises otherwise
  asl->i.c_vars_ = asl->i.o_vars_ = asl->i.n_var_;

  // TODO: allocate arrays as fg_read does
  int nv1 = asl->i.n_var_ + asl->i.nsufext[ASL_Sufkind_var];
  int ncom = 0;
  int nv = nv1 + ncom;
  int nc0 = asl->i.n_con_;
  int nc = nc0 + asl->i.nsufext[ASL_Sufkind_con];
  int no = asl->i.n_obj_;
  //int nvc = asl->i.c_vars_;
  //int nvo = asl->i.o_vars_;
  int nlcon = asl->i.n_lcon_;
  int nco = nc + no + nlcon;
  asl->i.ncom0_ = asl->i.combc_ + asl->i.como_;
  asl->i.ncom1_ = asl->i.comc1_ + asl->i.como1_;
  unsigned x =
      nco * sizeof(cde) + no * sizeof(ograd*)
    + nv * (sizeof(expr_v) + 2 * sizeof(int))
    //+ asl->i.ncom0_ * sizeof(cexp)
    + asl->i.ncom1_ * sizeof(cexp1)
    //+ nfunc * sizeof(func_info*)
    //+ nvref * sizeof(int)
    + no;
  expr_v *e = asl->I.var_e_ =
      reinterpret_cast<expr_v*>(M1zapalloc_ASL(&asl_->i, x));
  for (int i = 0; i < h.num_vars; ++i, ++e) {
    e->op = r_ops_[OPVARVAL];
    e->a = i;
  }
}

ExprFactory::~ExprFactory() {
  ASL_free(&asl_);
}

NumericConstant ExprFactory::CreateNumericConstant(double value) {
  expr_n *e = reinterpret_cast<expr_n*>(mem_ASL(asl_, asl_->i.size_expr_n_));
  e->op = reinterpret_cast<efunc_n*>(r_ops_[OPNUM]);
  e->v = value;
  return Expr::Create<NumericConstant>(reinterpret_cast<expr*>(e));
}

Variable ExprFactory::CreateVariable(int var_index) {
  assert(var_index >= 0 && var_index < asl_->i.n_var_);
  return Expr::Create<Variable>(reinterpret_cast<expr*>(
      reinterpret_cast<ASL_fg*>(asl_)->I.var_e_ + var_index));
}
}

