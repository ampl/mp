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

#include "solvers/util/nl.h"

#include <algorithm>
#include <cstring>

extern "C" void bswap_ASL(void *x, unsigned long L);

namespace {
static double DVALUE[] = {
#include "dvalue.hd"
};

// Initialize count variables in the array v starting from offset by
// setting their variable indices to var_index.
void InitVars(expr_v *v, int offset, int count, int var_index) {
  offset -= count;
  while (--count >= 0)
    v[offset + count].a = var_index;
}
}

namespace ampl {

internal::ASLBuilder::ASLBuilder(ASL &asl, const char *stub, const NLHeader &h)
: asl_(asl), r_ops_(0), nv1_(0), nz_(0), nderp_(0) {
  std::size_t stub_len = std::strlen(stub);
  Edaginfo &info = asl.i;
  info.filename_ = reinterpret_cast<char*>(M1alloc_ASL(&info, stub_len + 5));
  std::strcpy(info.filename_, stub);
  info.stub_end_ = info.filename_ + stub_len;
  std::strcpy(info.filename_ + stub_len, ".nl");

  info.binary_nl_ = h.format;
  if (h.format == NLHeader::BINARY_SWAPPED) {
    info.binary_nl_ = 1 << (3 - Arith_Kind_ASL);
    info.iadjfcn = info.dadjfcn = bswap_ASL;
  }
  info.xscanf_ = info.binary_nl_ ? bscanf : ascanf;

  info.ampl_options_[0] = h.num_options;
  for (int i = 0; i < ampl::MAX_NL_OPTIONS; ++i)
    info.ampl_options_[i + 1] = h.options[i];
  info.ampl_vbtol_ = h.ampl_vbtol;

  info.n_var_ = h.num_vars;
  info.n_con_ = h.num_algebraic_cons;
  info.n_obj_ = h.num_objs;
  info.nranges_ = h.num_ranges;
  info.n_eqn_ = h.num_eqns;
  info.n_lcon_ = h.num_logical_cons;

  info.nlc_ = h.num_nl_cons;
  info.nlo_ = h.num_nl_objs;
  info.n_cc_ = h.num_compl_conds;
  info.nlcc_ = h.num_nl_compl_conds;
  info.ndcc_ = h.num_compl_dbl_ineqs;
  info.nzlb_ = h.num_compl_vars_with_nz_lb;

  info.nlnc_ = h.num_nl_net_cons;
  info.lnc_ = h.num_linear_net_cons;

  info.nlvc_ = h.num_nl_vars_in_cons;
  info.nlvo_ = h.num_nl_vars_in_objs;
  info.nlvb_ = h.num_nl_vars_in_both;

  info.nwv_ = h.num_linear_net_vars;
  info.nfunc_ = h.num_funcs;
  info.flags = h.flags;

  info.nbv_ = h.num_linear_binary_vars;
  info.niv_ = h.num_linear_integer_vars;
  info.nlvbi_ = h.num_nl_integer_vars_in_both;
  info.nlvci_ = h.num_nl_integer_vars_in_cons;
  info.nlvoi_ = h.num_nl_integer_vars_in_objs;

  info.nzc_ = h.num_con_nonzeros;
  info.nzo_ = h.num_obj_nonzeros;

  info.maxrownamelen_ = h.max_con_name_len;
  info.maxcolnamelen_ = h.max_var_name_len;

  info.comb_ = h.num_common_exprs_in_both;
  info.comc_ = h.num_common_exprs_in_cons;
  info.como_ = h.num_common_exprs_in_objs;
  info.comc1_ = h.num_common_exprs_in_cons1;
  info.como1_ = h.num_common_exprs_in_objs1;

  info.nclcon_ = info.n_con_ + info.n_lcon_;

  if (h.num_algebraic_cons < 0 || h.num_vars <= 0 || h.num_objs < 0) {
    throw ASLError(ASL_readerr_corrupt,
        fmt::format("invalid problem dimensions: M = {}, N = {}, NO = {}",
          h.num_algebraic_cons, h.num_vars, h.num_objs));
  }

  info.n_var0 = info.n_var1 = info.n_var_;
  info.n_con0 = info.n_con1 = info.n_con_;
  info.x0len_ = std::max(std::max(info.nlvo_, info.nlvc_), 1) * sizeof(double);
  info.x0kind_ = ASL_first_x;
  info.n_conjac_[0] = 0;
  info.n_conjac_[1] = info.n_con_;

  info.c_vars_ = info.o_vars_ = info.n_var_;
}

void internal::ASLBuilder::BeginBuild(int flags) {
  bool linear = asl_.i.ASLtype == ASL_read_f;

  // Includes allocation of LUv, LUrhs, A_vals or Cgrad, etc.
  flagsave_ASL(&asl_, flags);

  Edaginfo &info = asl_.i;
  int nlcon = info.n_lcon_;
  if (nlcon && (flags & ASL_allow_CLP) == 0)
    throw ASLError(ASL_readerr_CLP, "cannot handle logical constraints");

  // TODO: test ASLBuilder
  bool readall = (flags & ASL_keep_all_suffixes) != 0;
  if (!info.size_expr_n_)
    info.size_expr_n_ = sizeof(expr_n);

  Edag1info &info1 = reinterpret_cast<ASL_fg&>(asl_).I;
  int ncom = 0;
  if (!linear) {
    r_ops_ = info1.r_ops_;
    if (!r_ops_)
      r_ops_ = r_ops_ASL;
    if (info.c_cexp1st_)
      *info.c_cexp1st_ = 0;
    if (info.o_cexp1st_)
      *info.o_cexp1st_ = info.comc1_;
    if (info.nfunc_)
      func_add(&asl_);
    ncom = info.comb_ + info.comc_ + info.como_ + info.comc1_ + info.como1_;
  }

  int nc0 = info.n_con_;
  int nc = nc0 + info.nsufext[ASL_Sufkind_con];
  int no = info.n_obj_;
  int nvc = info.c_vars_;
  int nvo = info.o_vars_;
  int nco = nc + no + nlcon;
  if (no < 0 || nco <= 0) {
    throw ASLError(ASL_readerr_corrupt,
        fmt::format("invalid problem dimensions: "
            "nc = {}, no = {}, nlcon = {}", nc0, no, nlcon));
  }
  if (info.pi0_) {
    std::memset(info.pi0_, 0, nc * sizeof(double));
    if (info.havepi0_)
      std::memset(info.havepi0_, 0, nc);
  }
  int nxv = info.nsufext[ASL_Sufkind_var];
  int nvr = info.n_var_;  // nv for reading
  int nv0 = nvr + nxv;
  nv1_ = nv0;
  int nv = nv1_;
  unsigned x = 0, nvref = 0, maxfwd1 = 0;
  if (linear) {
    x = nco * sizeof(cde) + no * sizeof(ograd*) + nv * sizeof(expr_v) + no;
  } else {
    int max_var = nv = nv1_ + ncom;
    info.combc_ = info.comb_ + info.comc_;
    int ncom_togo = info.ncom0_ = info.combc_ + info.como_;
    int nzclim = info.ncom0_ >> 3;
    info.ncom1_ = info.comc1_ + info.como1_;
    int nv0b = nv1_ + info.comb_;
    int nv0c = nv0b + info.comc_;
    int nv01 = nv1_ + info.ncom0_;
    int nv011 = nv01 - 1;
    int last_cex = nv011;
    int lasta00 = nv1_ + 1;
    int lasta = lasta00, lasta0 = lasta00;
    info.amax_ = lasta;
    maxfwd1 = asl_.p.maxfwd_ + 1;
    nvref = 0;
    if (maxfwd1 > 1) {
      nvref = maxfwd1 * ((info.ncom0_ < asl_.p.vrefGulp_ ?
          info.ncom0_ : asl_.p.vrefGulp_) + 1);
    }
    x = nco * sizeof(cde) + no * sizeof(ograd*)
        + nv * (sizeof(expr_v) + 2 * sizeof(int)) + info.ncom0_ * sizeof(cexp)
        + info.ncom1_ * sizeof(cexp1) + info.nfunc_ * sizeof(func_info*)
        + nvref * sizeof(int) + no;
    int nvar0 = info.n_var0;
    int nvinc = info.n_var_ - nvar0 + nxv;
    if (!nvinc)
      nvar0 += info.ncom0_ + info.ncom1_;
  }

  if (info.X0_)
    std::memset(info.X0_, 0, nv1_ * sizeof(double));
  if (info.havex0_)
    std::memset(info.havex0_, 0, nv1_);
  expr_v *e = info1.var_e_ =
      reinterpret_cast<expr_v*>(M1zapalloc_ASL(&info, x));
  info1.con_de_ = reinterpret_cast<cde*>(e + nv);
  info1.lcon_de_ = info1.con_de_ + nc;
  info1.obj_de_ = info1.lcon_de_ + nlcon;
  info.Ograd_ = reinterpret_cast<ograd**>(info1.obj_de_ + no);

  if (linear) {
    info.objtype_ = reinterpret_cast<char*>(info.Ograd_ + no);
  } else {
    info1.var_ex_ = e + nv1_;
    info1.var_ex1_ = info1.var_ex_ + info.ncom0_;
    for (int k = 0; k < nv; ++e, ++k) {
      e->op = r_ops_[OPVARVAL];
      e->a = k;
    }
    if (info.skip_int_derivs_) {
      if (info.nlvbi_)
        InitVars(e, info.nlvb_, info.nlvbi_, nv1_);
      if (info.nlvci_)
        InitVars(e, info.nlvb_ + info.nlvc_, info.nlvci_, nv1_);
      if (info.nlvoi_)
        InitVars(e, info.nlvb_ + info.nlvc_ + info.nlvo_, info.nlvoi_, nv1_);
    }
    info1.cexps_ = reinterpret_cast<cexp*>(info.Ograd_ + no);
    info1.cexps1_ = reinterpret_cast<cexp1*>(info1.cexps_ + info.ncom0_);
    info.funcs_ = reinterpret_cast<func_info**>(info1.cexps1_
        + info.ncom1_);
    int *zc = reinterpret_cast<int*>(info.funcs_ + info.nfunc_);
    int *zci = zc + nv;
    int *vrefx = zci + nv;
    info.objtype_ = reinterpret_cast<char*>(vrefx + nvref);
    if (nvref) {
      int *vrefnext = vrefx + maxfwd1;
      nvref -= maxfwd1;
    }
    int last_d = 0;
  }
  if (info.n_cc_ && !info.cvar_)
    info.cvar_ = reinterpret_cast<int*>(M1alloc_ASL(&info, nc * sizeof(int)));
  if (info.cvar_)
    std::memset(info.cvar_, 0, nc * sizeof(int));
  int *ka = 0;
  nz_ = 0;
  nderp_ = 0;
}

void internal::ASLBuilder::EndBuild() {
  bool linear = asl_.i.ASLtype == ASL_read_f;
  Edaginfo &info = asl_.i;
  if (!linear) {
    // Make amax long enough for nlc to handle
    // var_e[i].a for common variables i.
    if (info.ncom0_) {
      int i = info.comb_ + info.como_;
      if (i < info.combc_)
        i = info.combc_;
      if ((i += nv1_ + 1) > info.amax_)
        info.amax_ = i;
    }
    info.adjoints_ = reinterpret_cast<double*>(
        M1zapalloc_ASL(&asl_.i, info.amax_ * Sizeof(double)));
    info.adjoints_nv1_ = &info.adjoints_[nv1_];
    info.nderps_ += nderp_;
  }
  // TODO
  //adjust(S, flags);
  info.nzjac_ = nz_;
  if (!info.Lastx_) {
    info.Lastx_ = reinterpret_cast<double*>(
        M1alloc_ASL(&info, nv1_ * sizeof(double)));
  }
  if (!linear) {
    asl_.p.Objval = asl_.p.Objval_nomap = obj1val_ASL;
    asl_.p.Objgrd = asl_.p.Objgrd_nomap = obj1grd_ASL;
    asl_.p.Conval = con1val_ASL;
    asl_.p.Jacval = jac1val_ASL;
    asl_.p.Conival = asl_.p.Conival_nomap = con1ival_ASL;
    asl_.p.Congrd = asl_.p.Congrd_nomap = con1grd_ASL;
    asl_.p.Lconval = lcon1val_ASL;
    asl_.p.Xknown = x1known_ASL;
  }
  prob_adj_ASL(&asl_);
}

ExprFactory::ExprFactory(const NLHeader &h, const char *stub, int flags)
: asl_(ASL_alloc(ASL_read_fg)) {
  internal::ASLBuilder builder(*asl_, stub, h);

  // TODO: this part should be optional
  ASL_fg *asl = reinterpret_cast<ASL_fg*>(asl_);
  for (int i = 0; i < N_OPS; ++i)
    r_ops_[i] = reinterpret_cast<efunc*>(i);
  asl->I.r_ops_ = r_ops_;

  builder.BeginBuild(flags);
  builder.EndBuild();
}

ExprFactory::~ExprFactory() {
  ASL_free(&asl_);
}

template <typename ExprT>
ExprT ExprFactory::MakeExpr(int opcode, NumericExpr lhs, NumericExpr rhs) {
  expr *e = Allocate<expr>();
  e->op = reinterpret_cast<efunc*>(r_ops_[opcode]);
  e->L.e = lhs.expr_;
  e->R.e = rhs.expr_;
  e->a = asl_->i.n_var_ + asl_->i.nsufext[ASL_Sufkind_var];
  e->dL = DVALUE[opcode];  // for UMINUS, FLOOR, CEIL
  return Expr::Create<ExprT>(e);
}

UnaryExpr ExprFactory::MakeUnaryExpr(int opcode, NumericExpr arg) {
  CheckOpCode(opcode, Expr::UNARY, "unary");
  UnaryExpr expr = MakeExpr<UnaryExpr>(opcode, arg, NumericExpr());
  expr.expr_->dL = DVALUE[opcode];  // for UMINUS, FLOOR, CEIL
  return expr;
}

BinaryExpr ExprFactory::MakeBinaryExpr(
    int opcode, NumericExpr lhs, NumericExpr rhs) {
  CheckOpCode(opcode, Expr::BINARY, "binary");
  BinaryExpr expr = MakeExpr<BinaryExpr>(opcode, lhs, rhs);
  expr.expr_->dL = 1;
  expr.expr_->dR = DVALUE[opcode];  // for PLUS, MINUS, REM
  return expr;
}

VarArgExpr ExprFactory::MakeVarArgExpr(
    int opcode, int num_args, NumericExpr *args) {
  assert(num_args >= 0);
  CheckOpCode(opcode, Expr::VARARG, "vararg");
  expr_va *result = Allocate<expr_va>();
  result->op = reinterpret_cast<efunc *>(r_ops_[opcode]);
  de *d = result->L.d = Allocate<de>(num_args * sizeof(de) + sizeof(expr*));
  for (int i = 0; i < num_args; ++i)
    d[i].e = args[i].expr_;
  d[num_args].e = 0;
  return Expr::Create<VarArgExpr>(reinterpret_cast<expr*>(result));
}

Variable ExprFactory::MakeVariable(int var_index) {
  assert(var_index >= 0 && var_index < asl_->i.n_var_);
  return Expr::Create<Variable>(
      reinterpret_cast<expr*>(reinterpret_cast<ASL_fg*>(asl_)->I.var_e_
          + var_index));
}
}
