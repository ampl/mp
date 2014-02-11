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

#include <cstring>

extern "C" void bswap_ASL(void *x, unsigned long L);

namespace {
// Initialize count variables in the array v starting from offset by
// setting their variable indices to var_index.
void InitVars(expr_v *v, int offset, int count, int var_index) {
  offset -= count;
  while (--count >= 0)
    v[offset + count].a = var_index;
}
}

namespace ampl {

void internal::InitASL(ASL &asl, const char *stub, const NLHeader &h) {
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

  info.n_var0 = info.n_var1 = info.n_var_;
  info.n_con0 = info.n_con1 = info.n_con_;
  int nlv = info.nlvc_;
  if (nlv < info.nlvo_)
    nlv = info.nlvo_;
  if (nlv <= 0)
    nlv = 1;
  info.x0len_ = nlv * sizeof(double);
  info.x0kind_ = ASL_first_x;
  info.n_conjac_[0] = 0;
  info.n_conjac_[1] = info.n_con_;

  info.c_vars_ = info.o_vars_ = info.n_var_;
}

class ReadError: public ampl::Error {
private:
  int error_code_;

public:
  ReadError(int error_code, fmt::StringRef message) :
      Error(message), error_code_(error_code) {
  }
};

ExprFactory::ExprFactory(const NLHeader &h, const char *stub) :
    asl_(ASL_alloc(ASL_read_fg)) {
  internal::InitASL(*asl_, stub, h);

  int flags = 0; // TODO: flags (ASL_return_read_err, ASL_allow_CLP)

  // Includes allocation of LUv, LUrhs, A_vals or Cgrad, etc.
  flagsave_ASL(asl_, flags);
  int nlcon = asl_->i.n_lcon_;
  if (nlcon && (flags & ASL_allow_CLP) == 0)
    throw ReadError(ASL_readerr_CLP, "cannot handle logical constraints");

  bool readall = (flags & ASL_keep_all_suffixes) != 0;
  if (!asl_->i.size_expr_n_)
    asl_->i.size_expr_n_ = sizeof(expr_n);

  asl_->i.xscanf_ = asl_->i.binary_nl_ ? bscanf : ascanf;

  bool just_linear = false;
  ASL_fg *asl = reinterpret_cast<ASL_fg*>(asl_);
  int ncom = 0;
  efunc **r_ops = 0;
  if (!just_linear) {
    r_ops = asl->I.r_ops_;
    if (!r_ops)
      r_ops = r_ops_ASL;
    if (asl_->i.c_cexp1st_)
      *asl_->i.c_cexp1st_ = 0;
    if (asl_->i.o_cexp1st_)
      *asl_->i.o_cexp1st_ = asl_->i.comc1_;
    if (asl_->i.nfunc_)
      func_add(asl_);
    ncom = asl_->i.comb_ + asl_->i.comc_ + asl_->i.como_ + asl_->i.comc1_
        + asl_->i.como1_;
  }

  int nc0 = asl_->i.n_con_;
  int nc = nc0 + asl_->i.nsufext[ASL_Sufkind_con];
  int no = asl_->i.n_obj_;
  int nvc = asl_->i.c_vars_;
  int nvo = asl_->i.o_vars_;
  int nco = nc + no + nlcon;
  // TODO
  //if (no < 0 || nco <= 0)
  //  throw ReadError(ASL_readerr_corrupt, "corrupt .nl file");
  if (asl_->i.pi0_) {
    memset(asl_->i.pi0_, 0, nc * sizeof(double));
    if (asl_->i.havepi0_)
      memset(asl_->i.havepi0_, 0, nc);
  }
  int nxv = asl_->i.nsufext[ASL_Sufkind_var];
  int nvr = asl_->i.n_var_;  // nv for reading
  int nv0 = nvr + nxv;
  int nv1 = nv0;
  int nv = nv1;
  unsigned x = 0, nvref = 0, maxfwd1 = 0;
  if (just_linear) {
    x = nco * sizeof(cde) + no * sizeof(ograd*) + nv * sizeof(expr_v) + no;
  } else {
    int max_var = nv = nv1 + ncom;
    asl->i.combc_ = asl->i.comb_ + asl->i.comc_;
    int ncom_togo = asl->i.ncom0_ = asl->i.combc_ + asl->i.como_;
    int nzclim = asl->i.ncom0_ >> 3;
    asl->i.ncom1_ = asl->i.comc1_ + asl->i.como1_;
    int nv0b = nv1 + asl->i.comb_;
    int nv0c = nv0b + asl->i.comc_;
    int nv01 = nv1 + asl->i.ncom0_;
    int nv011 = nv01 - 1;
    int last_cex = nv011;
    int lasta00 = nv1 + 1;
    int lasta = lasta00, lasta0 = lasta00;
    asl->i.amax_ = lasta;
    maxfwd1 = maxfwd + 1;
    nvref = 0;
    if (maxfwd1 > 1) {
      nvref = maxfwd1
          * ((asl->i.ncom0_ < vrefGulp ? asl->i.ncom0_ : vrefGulp) + 1);
    }
    x = nco * sizeof(cde) + no * sizeof(ograd*)
        + nv * (sizeof(expr_v) + 2 * sizeof(int)) + asl->i.ncom0_ * sizeof(cexp)
        + asl->i.ncom1_ * sizeof(cexp1) + asl->i.nfunc_ * sizeof(func_info*)
        + nvref * sizeof(int) + no;
    int nvar0 = asl_->i.n_var0;
    int nvinc = asl_->i.n_var_ - nvar0 + nxv;
    if (!nvinc)
      nvar0 += asl_->i.ncom0_ + asl_->i.ncom1_;
  }

  if (asl_->i.X0_)
    memset(asl_->i.X0_, 0, nv1 * sizeof(double));
  if (asl_->i.havex0_)
    memset(asl_->i.havex0_, 0, nv1);
  expr_v *e = asl->I.var_e_ = reinterpret_cast<expr_v*>(M1zapalloc(x));
  asl->I.con_de_ = reinterpret_cast<cde*>(e + nv);
  asl->I.lcon_de_ = asl->I.con_de_ + nc;
  asl->I.obj_de_ = asl->I.lcon_de_ + nlcon;
  asl->i.Ograd_ = reinterpret_cast<ograd**>(asl->I.obj_de_ + no);

  if (just_linear) {
    asl_->i.objtype_ = reinterpret_cast<char*>(asl_->i.Ograd_ + no);
  } else {
    var_ex = e + nv1;
    var_ex1 = var_ex + asl_->i.ncom0_;
    for (int k = 0; k < nv; ++e, ++k) {
      e->op = r_ops_[OPVARVAL];
      e->a = k;
    }
    if (asl->i.skip_int_derivs_) {
      if (asl_->i.nlvbi_)
        InitVars(var_e, asl_->i.nlvb_, asl_->i.nlvbi_, nv1);
      if (asl_->i.nlvci_)
        InitVars(var_e, asl_->i.nlvb_ + asl_->i.nlvc_, asl_->i.nlvci_, nv1);
      if (asl_->i.nlvoi_) {
        InitVars(var_e, asl_->i.nlvb_ + asl_->i.nlvc_ + asl_->i.nlvo_,
            asl_->i.nlvoi_, nv1);
      }
    }
    asl->I.cexps_ = reinterpret_cast<cexp*>(asl_->i.Ograd_ + no);
    asl->I.cexps1_ = reinterpret_cast<cexp1*>(asl->I.cexps_ + asl_->i.ncom0_);
    asl_->i.funcs_ = reinterpret_cast<func_info**>(asl->I.cexps1_
        + asl_->i.ncom1_);
    int *zc = reinterpret_cast<int*>(asl_->i.funcs_ + asl_->i.nfunc_);
    int *zci = zc + nv;
    int *vrefx = zci + nv;
    asl_->i.objtype_ = reinterpret_cast<char*>(vrefx + nvref);
    if (nvref) {
      int *vrefnext = vrefx + maxfwd1;
      nvref -= maxfwd1;
    }
    int last_d = 0;
  }
  if (asl_->i.n_cc_ && !asl_->i.cvar_)
    asl_->i.cvar_ = reinterpret_cast<int*>(M1alloc(nc * sizeof(int)));
  if (asl_->i.cvar_)
    std::memset(asl_->i.cvar_, 0, nc * sizeof(int));
  int *ka = 0;
  int nz = 0;
  int nderp = 0;

  // TODO: this part should be optional
  for (int i = 0; i < N_OPS; ++i)
    r_ops_[i] = reinterpret_cast<efunc*>(i);
  asl->I.r_ops_ = r_ops_;

  // TODO: this part should be executed after problem construction
  if (!just_linear) {
    // Make amax long enough for nlc to handle
    // var_e[i].a for common variables i.
    if (asl_->i.ncom0_) {
      int i = asl_->i.comb_ + asl_->i.como_;
      if (i < asl_->i.combc_)
        i = asl_->i.combc_;
      if ((i += nv1 + 1) > asl_->i.amax_)
        asl_->i.amax_ = i;
    }
    asl_->i.adjoints_ = reinterpret_cast<double*>(
        M1zapalloc(asl_->i.amax_ * Sizeof(double)));
    asl_->i.adjoints_nv1_ = &asl_->i.adjoints_[nv1];
    asl_->i.nderps_ += nderp;
  }
  // TODO
  if (!just_linear) {
    asl_->p.Objval = asl_->p.Objval_nomap = obj1val_ASL;
    asl_->p.Objgrd = asl_->p.Objgrd_nomap = obj1grd_ASL;
    asl_->p.Conval = con1val_ASL;
    asl_->p.Jacval = jac1val_ASL;
    asl_->p.Conival = asl_->p.Conival_nomap = con1ival_ASL;
    asl_->p.Congrd = asl_->p.Congrd_nomap = con1grd_ASL;
    asl_->p.Lconval = lcon1val_ASL;
    asl_->p.Xknown = x1known_ASL;
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
  return Expr::Create<Variable>(
      reinterpret_cast<expr*>(reinterpret_cast<ASL_fg*>(asl_)->I.var_e_
          + var_index));
}
}

