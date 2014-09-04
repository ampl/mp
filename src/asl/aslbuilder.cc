/*
 An ASL problem builder.

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

#define ASL_PRESERVE_DEFINES
#include "aslbuilder.h"

#include "mp/nl.h"
#include "opcode.hd"

#include <algorithm>
#include <cstring>

// Include fg_read.c to access static functions defined there.
#define fg_read_ASL fg_read2_ASL
#include "fg_read.c"
#undef asl
#undef max

extern "C" void bswap_ASL(void *x, unsigned long L);

namespace {

// If Double_Align is defined, bring size up to the nearest
// multiple of sizeof(long).
inline int AddPadding(int size) {
#ifdef Double_Align
  return size;
#else
  return (size + sizeof(long) - 1) & ~(sizeof(long) - 1);
#endif
}

// Initialize count variables in the array v starting from offset by
// setting their variable indices to var_index.
void InitVars(expr_v *v, int offset, int count, int var_index) {
  offset -= count;
  while (--count >= 0)
    v[offset + count].a = var_index;
}

double MissingFunc(arglist *al) {
  static char error[] = "attempt to call unavailable function";
  al->Errmsg = error;
  return 0;
}
}  // namespace

namespace mp {
namespace internal {

const double ASLBuilder::DVALUE[] = {
#include "dvalue.hd"
};

template <typename T>
inline T *ASLBuilder::ZapAllocate(std::size_t size) {
  return reinterpret_cast<T*>(M1zapalloc_ASL(&asl_->i, size));
}

template <typename T>
T *ASLBuilder::AllocateSuffixValues(
    T *&values, int num_values, int nx, int nx1) {
  if (!values)
    values = Allocate<T>(nx1 * sizeof(T));
  if (num_values < nx)
    std::memset(values, 0, nx * sizeof(T));
  if (nx < nx1)
    std::memset(values + nx, 0, (nx1 - nx) * sizeof(T));
  return values;
}

void ASLBuilder::SetObjOrCon(
    int index, cde *d, int *cexp1_end, ::expr *e, int **z) {
  bool linear = asl_->i.ASLtype == ASL_read_f;
  ASL_fg *asl = reinterpret_cast<ASL_fg*>(asl_);
  if (!linear) {
    static_->_lastc1 = static_->_last_cex - static_->_nv011;
    if (cexp1_end)
      cexp1_end[index + 1] = static_->_lastc1;
    if (static_->_amax1 < static_->_lasta)
      static_->_amax1 = static_->_lasta;
    if (static_->_co_first) {
      static_->_co_first = 0;
      if (static_->_imap_len < static_->_lasta)
        imap_alloc(static_);
      asl->I.f_b_ = funnelfix(asl->I.f_b_);
      asl->I.f_c_ = funnelfix(asl->I.f_c_);
      asl->I.f_o_ = funnelfix(asl->I.f_o_);
    }
    if (!static_->_lastj) {
      static_->_lasta = static_->_lasta0;
      static_->_last_d = 0;
    }
    static_->_lastj = 0;
  }
  d += index;
  d->e = e;
  if (!linear) {
    int alen = static_->_lasta - static_->_lasta0;
    if (static_->_imap_len < static_->_lasta)
      imap_alloc(static_);
    if (z) {
      z += index;
      *z = 0;
    }
    comsubs(static_, alen, d, z);
    static_->_firstc1 = static_->_lastc1;
  }
}

::expr *ASLBuilder::MakeConstant(double value) {
  expr_n *result = Allocate<expr_n>(asl_->i.size_expr_n_);
  result->op = reinterpret_cast<efunc_n*>(OPNUM);
  result->v = value;
  return reinterpret_cast< ::expr*>(result);
}

::expr *ASLBuilder::DoMakeUnary(expr::Kind kind, Expr arg) {
  ::expr *e = Allocate< ::expr>();
  int opcode = expr::opcode(kind);
  e->op = reinterpret_cast<efunc*>(opcode);
  e->L.e = arg.expr_;
  e->a = asl_->i.n_var_ + asl_->i.nsufext[ASL_Sufkind_var];
  e->dL = DVALUE[opcode];  // for UMINUS, FLOOR, CEIL
  return e;
}

::expr *ASLBuilder::DoMakeBinary(expr::Kind kind, Expr lhs, Expr rhs) {
  ::expr *e = Allocate< ::expr>();
  int opcode = expr::opcode(kind);
  e->op = reinterpret_cast<efunc*>(opcode);
  e->L.e = lhs.expr_;
  e->R.e = rhs.expr_;
  e->a = asl_->i.n_var_ + asl_->i.nsufext[ASL_Sufkind_var];
  e->dL = 1;
  e->dR = DVALUE[opcode];  // for PLUS, MINUS, REM
  return e;
}

::expr *ASLBuilder::MakeIf(
    expr::Kind kind, LogicalExpr condition, Expr true_expr, Expr false_expr) {
  expr_if *e = Allocate<expr_if>();
  e->op = r_ops_[opcode(kind)];
  e->e = condition.expr_;
  e->T = true_expr.expr_;
  e->F = false_expr.expr_;
  return reinterpret_cast< ::expr*>(e);
}

::expr *ASLBuilder::MakeIterated(expr::Kind kind, ArrayRef<Expr> args) {
  int num_args = SafeInt<int>(args.size()).value();
  ::expr *e = Allocate< ::expr>(
      sizeof(::expr) - sizeof(double) +
      SafeInt<int>(num_args + 1) * sizeof(::expr*));
  e->op = reinterpret_cast<efunc*>(opcode(kind));
  e->L.ep = reinterpret_cast< ::expr**>(&e->dR);
  e->R.ep = e->L.ep + num_args;
  ::expr **arg_ptrs = e->L.ep;
  for (int i = 0; i < num_args; ++i)
    arg_ptrs[i] = args[i].expr_;
  return e;
}

ASLBuilder::ASLBuilder(ASL *asl)
: asl_(asl), own_asl_(false), r_ops_(0), flags_(0),
  nz_(0), nderp_(0), static_(0) {
  if (!asl) {
    asl_ = ASL_alloc(ASL_read_fg);
    own_asl_ = true;
  }
}

ASLBuilder::~ASLBuilder() {
  if (own_asl_)
    ASL_free(&asl_);
  if (static_)
    delete static_;
}

void ASLBuilder::set_stub(const char *stub) {
  std::size_t stub_len = std::strlen(stub);
  Edaginfo &info = asl_->i;
  info.filename_ = reinterpret_cast<char*>(M1alloc_ASL(&info, stub_len + 5));
  std::strcpy(info.filename_, stub);
  info.stub_end_ = info.filename_ + stub_len;
  std::strcpy(info.filename_ + stub_len, ".nl");
}

void ASLBuilder::InitASL(const NLHeader &h) {
  if (!static_)
    static_ = new Static();
  Edaginfo &info = asl_->i;
  info.binary_nl_ = h.format;
  if (h.arith_kind != arith::UNKNOWN) {
    arith::Kind arith_kind = arith::GetKind();
    if (arith_kind != h.arith_kind &&
        arith::IsIEEE(arith_kind) && arith::IsIEEE(h.arith_kind)) {
      info.binary_nl_ = h.arith_kind << 1;
      info.iadjfcn = info.dadjfcn = bswap_ASL;
    }
  }
  info.xscanf_ = info.binary_nl_ ? bscanf : ascanf;

  info.ampl_options_[0] = h.num_options;
  for (int i = 0; i < mp::MAX_NL_OPTIONS; ++i)
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
  info.comc1_ = h.num_common_exprs_in_single_cons;
  info.como1_ = h.num_common_exprs_in_single_objs;

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

void ASLBuilder::BeginBuild(const NLHeader &h) {
  InitASL(h);

  bool linear = asl_->i.ASLtype == ASL_read_f;

  Edaginfo &info = asl_->i;
  if ((flags_ & ASL_COLUMNWISE) != 0) {
    info.A_vals_ = reinterpret_cast<double*>(
        Malloc(info.nzc_ * sizeof(*info.A_vals_)));
    assert(asl_->i.A_vals_ != 0);
  }

  // Includes allocation of LUv, LUrhs, A_vals or Cgrad, etc.
  flagsave_ASL(asl_, flags_);
  ed_reset(static_, asl_);

  int nlcon = info.n_lcon_;
  if (nlcon && (flags_ & ASL_allow_CLP) == 0)
    throw ASLError(ASL_readerr_CLP, "cannot handle logical constraints");

  // TODO: test
  bool readall = (flags_ & ASL_keep_all_suffixes) != 0;
  if (!info.size_expr_n_)
    info.size_expr_n_ = sizeof(expr_n);

  Edag1info &info1 = reinterpret_cast<ASL_fg*>(asl_)->I;
  int ncom = 0;
  if (!linear) {
    if ((flags_ & ASL_STANDARD_OPCODES) != 0) {
      for (int i = 0; i < N_OPS; ++i)
        standard_opcodes_[i] = reinterpret_cast<efunc*>(i);
      r_ops_ = standard_opcodes_;
    } else {
      r_ops_ = info1.r_ops_;
      if (!r_ops_)
        r_ops_ = r_ops_ASL;
    }
    if (info.c_cexp1st_)
      *info.c_cexp1st_ = 0;
    if (info.o_cexp1st_)
      *info.o_cexp1st_ = info.comc1_;
    if (info.nfunc_)
      func_add(asl_);
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
  static_->_nv0 = nvr + nxv;
  static_->_nv1 = static_->_nv0;
  int nv = static_->_nv1;
  unsigned x = 0, maxfwd1 = 0;
  static_->_nvref = 0;
  if (linear) {
    x = nco * sizeof(cde) + no * sizeof(ograd*) + nv * sizeof(expr_v) + no;
  } else {
    static_->_max_var = nv = static_->_nv1 + ncom;
    info.combc_ = info.comb_ + info.comc_;
    static_->_ncom_togo = info.ncom0_ = info.combc_ + info.como_;
    static_->_nzclim = info.ncom0_ >> 3;
    info.ncom1_ = info.comc1_ + info.como1_;
    static_->_nv0b = static_->_nv1 + info.comb_;
    static_->_nv0c = static_->_nv0b + info.comc_;
    static_->_nv01 = static_->_nv1 + info.ncom0_;
    static_->_nv011 = static_->_nv01 - 1;
    static_->_last_cex = static_->_nv011;
    static_->_lasta00 = static_->_nv1 + 1;
    static_->_lasta = static_->_lasta00;
    static_->_lasta0 = static_->_lasta00;
    info.amax_ = static_->_lasta;
    maxfwd1 = asl_->p.maxfwd_ + 1;
    static_->_nvref = 0;
    if (maxfwd1 > 1) {
      static_->_nvref = maxfwd1 * ((info.ncom0_ < asl_->p.vrefGulp_ ?
          info.ncom0_ : asl_->p.vrefGulp_) + 1);
    }
    x = nco * sizeof(cde) + no * sizeof(ograd*)
        + nv * (sizeof(expr_v) + 2 * sizeof(int)) + info.ncom0_ * sizeof(cexp)
        + info.ncom1_ * sizeof(cexp1) + info.nfunc_ * sizeof(func_info*)
        + static_->_nvref * sizeof(int) + no;
    int nvar0 = info.n_var0;
    int nvinc = info.n_var_ - nvar0 + nxv;
    if (!nvinc)
      nvar0 += info.ncom0_ + info.ncom1_;
  }

  if (info.X0_)
    std::memset(info.X0_, 0, static_->_nv1 * sizeof(double));
  if (info.havex0_)
    std::memset(info.havex0_, 0, static_->_nv1);
  expr_v *e = info1.var_e_ = ZapAllocate<expr_v>(x);
  info1.con_de_ = reinterpret_cast<cde*>(e + nv);
  info1.lcon_de_ = info1.con_de_ + nc;
  info1.obj_de_ = info1.lcon_de_ + nlcon;
  info.Ograd_ = reinterpret_cast<ograd**>(info1.obj_de_ + no);

  if (linear) {
    info.objtype_ = reinterpret_cast<char*>(info.Ograd_ + no);
  } else {
    info1.var_ex_ = e + static_->_nv1;
    info1.var_ex1_ = info1.var_ex_ + info.ncom0_;
    for (int k = 0; k < nv; ++e, ++k) {
      e->op = r_ops_[OPVARVAL];
      e->a = k;
    }
    if (info.skip_int_derivs_) {
      if (info.nlvbi_)
        InitVars(e, info.nlvb_, info.nlvbi_, static_->_nv1);
      if (info.nlvci_)
        InitVars(e, info.nlvb_ + info.nlvc_, info.nlvci_, static_->_nv1);
      if (info.nlvoi_) {
        InitVars(e, info.nlvb_ + info.nlvc_ + info.nlvo_,
                 info.nlvoi_, static_->_nv1);
      }
    }
    info1.cexps_ = reinterpret_cast<cexp*>(info.Ograd_ + no);
    info1.cexps1_ = reinterpret_cast<cexp1*>(info1.cexps_ + info.ncom0_);
    info.funcs_ = reinterpret_cast<func_info**>(info1.cexps1_
        + info.ncom1_);
    static_->_zc = reinterpret_cast<int*>(info.funcs_ + info.nfunc_);
    static_->_zci = static_->_zc + nv;
    static_->_vrefx = static_->_zci + nv;
    info.objtype_ = reinterpret_cast<char*>(static_->_vrefx + static_->_nvref);
    if (static_->_nvref) {
      static_->_vrefnext = static_->_vrefx + maxfwd1;
      static_->_nvref -= maxfwd1;
    }
    static_->_last_d = 0;
  }
  if (info.n_cc_ && !info.cvar_)
    info.cvar_ = reinterpret_cast<int*>(M1alloc_ASL(&info, nc * sizeof(int)));
  if (info.cvar_)
    std::memset(info.cvar_, 0, nc * sizeof(int));
  int *ka = 0;
  nz_ = 0;
  nderp_ = 0;
}

void ASLBuilder::EndBuild() {
  bool linear = asl_->i.ASLtype == ASL_read_f;
  Edaginfo &info = asl_->i;
  if (!linear) {
    // Make amax long enough for nlc to handle
    // var_e[i].a for common variables i.
    if (info.ncom0_) {
      int i = info.comb_ + info.como_;
      if (i < info.combc_)
        i = info.combc_;
      if ((i += static_->_nv1 + 1) > info.amax_)
        info.amax_ = i;
    }
    info.adjoints_ = ZapAllocate<double>(info.amax_ * sizeof(double));
    info.adjoints_nv1_ = &info.adjoints_[static_->_nv1];
    info.nderps_ += nderp_;
  }
  // TODO
  //adjust(S, flags);
  info.nzjac_ = nz_;
  if (!info.Lastx_) {
    info.Lastx_ = reinterpret_cast<double*>(
        M1alloc_ASL(&info, static_->_nv1 * sizeof(double)));
  }
  if (!linear) {
    Edagpars &pars = asl_->p;
    pars.Objval = pars.Objval_nomap = obj1val_ASL;
    pars.Objgrd = pars.Objgrd_nomap = obj1grd_ASL;
    pars.Conval = con1val_ASL;
    pars.Jacval = jac1val_ASL;
    pars.Conival = pars.Conival_nomap = con1ival_ASL;
    pars.Congrd = pars.Congrd_nomap = con1grd_ASL;
    pars.Lconval = lcon1val_ASL;
    pars.Xknown = x1known_ASL;
  }
  prob_adj_ASL(asl_);
}

void ASLBuilder::SetObj(int index, obj::Type type, NumericExpr expr) {
  assert(0 <= index && index < asl_->i.n_obj_);
  asl_->i.objtype_[index] = type;
  SetObjOrCon(index, reinterpret_cast<ASL_fg*>(asl_)->I.obj_de_,
              asl_->i.o_cexp1st_, expr.expr_, asl_->i.zao_);
}

void ASLBuilder::SetCon(int index, NumericExpr expr) {
  assert(0 <= index && index < asl_->i.n_con_);
  SetObjOrCon(index, reinterpret_cast<ASL_fg*>(asl_)->I.con_de_,
              asl_->i.c_cexp1st_, expr.expr_, asl_->i.zac_);
}

void ASLBuilder::SetLogicalCon(int index, LogicalExpr expr) {
  assert(0 <= index && index < asl_->i.n_lcon_);
  SetObjOrCon(index, reinterpret_cast<ASL_fg*>(asl_)->I.lcon_de_,
              0, expr.expr_, 0);
}

Function ASLBuilder::AddFunction(
    const char *name, ufunc f, int num_args, func::Type type, void *info) {
  func_info *fi = func_lookup(asl_, name, 1);
  if (fi) {
    fi->funcp = f;
    fi->ftype = type;
    fi->nargs = num_args;
    fi->funcinfo = info;
    if (!asl_->i.funcsfirst_)
      asl_->i.funcsfirst_ = fi;
    else
      asl_->i.funcslast_->fnext = fi;
    asl_->i.funcslast_ = fi;
    fi->fnext = 0;
  }
  return Function(fi);
}

Function ASLBuilder::SetFunction(
    int index, fmt::StringRef name, int num_args, func::Type type) {
  assert(index >= 0 && index < asl_->i.nfunc_);
  func_info *fi = func_lookup_ASL(asl_, name.c_str(), 0);
  if (fi) {
    if (fi->nargs != num_args && fi->nargs >= 0 &&
        (num_args >= 0 || fi->nargs < -(num_args + 1))) {
      throw ASLError(ASL_readerr_argerr,
          fmt::format("function {}: disagreement of nargs: {} and {}",
                      name, fi->nargs, num_args));
    }
  } else {
    fi = Allocate<func_info>();
    fi->ftype = type;
    fi->nargs = num_args;
    fi->funcp = 0;
    int length = AddPadding(name.size() + 1);
    // TODO: name may not be null terminated
    fi->name = std::strcpy(Allocate<char>(length), name.c_str());
  }
  if (!fi->funcp && !(fi->funcp = dynlink(name.c_str()))) {
    if (!(flags_ & ASL_allow_missing_funcs)) {
      throw ASLError(ASL_readerr_unavail,
                     fmt::format("function {} not available", name));
    }
    fi->funcp = MissingFunc;
    fi->funcinfo = fi;
  }
  asl_->i.funcs_[index] = fi;
  return Function(fi);
}

ASLBuilder::SuffixHandler ASLBuilder::AddSuffix(
    int kind, int num_values, fmt::StringRef name) {
  bool readall = (flags_ & ASL_keep_all_suffixes) != 0;
  int item_type = kind & suf::MASK;
  SufDesc *d = 0;
  if (readall) {
    d = ZapAllocate<SufDesc>(sizeof(SufDesc) + name.size() + 1);
    d->next = asl_->i.suffixes[item_type];
    asl_->i.suffixes[item_type] = d;
    asl_->i.nsuff[item_type]++;
    asl_->i.nsuffixes++;
    std::copy(name.c_str(), name.c_str() + name.size(),
              d->sufname = reinterpret_cast<char*>(d + 1));
    d->kind = kind;
  } else {
    for (d = asl_->i.suffixes[item_type]; ; d = d->next) {
      if (!d)
        return SuffixHandler();  // Skip this suffix table.
      if (item_type == (d->kind & suf::MASK) &&
          !strncmp(name.c_str(), d->sufname, name.size()) &&
          !d->sufname[name.size()])
        if ((d->kind & suf::OUTONLY) != 0)
          return SuffixHandler();
        break;
    }
  }
  int nx = (&asl_->i.n_var_)[item_type];
  int nx1 = nx + d->nextra + asl_->i.nsufext[item_type];
  d->kind |= suf::INPUT;
  if ((d->kind & suf::FLOAT) != 0) {
    d->u.i = 0;
    return SuffixHandler(AllocateSuffixValues(d->u.r, num_values, nx, nx1));
  } else {
    d->u.r = 0;
    return SuffixHandler(AllocateSuffixValues(d->u.i, num_values, nx, nx1));
  }
  return SuffixHandler();
}

Variable ASLBuilder::MakeVariable(int var_index) {
  assert(var_index >= 0 && var_index < asl_->i.n_var_);
  return Expr::Create<Variable>(
      reinterpret_cast< ::expr*>(reinterpret_cast<ASL_fg*>(asl_)->I.var_e_
          + var_index));
}

UnaryExpr ASLBuilder::MakeUnary(expr::Kind kind, NumericExpr arg) {
  CheckKind<UnaryExpr>(kind, "unary");
  UnaryExpr expr = MakeUnary<UnaryExpr>(kind, arg);
  expr.expr_->dL = DVALUE[opcode(kind)];  // for UMINUS, FLOOR, CEIL
  return expr;
}

PiecewiseLinearExpr ASLBuilder::MakePiecewiseLinear(
    int num_breakpoints, const double *breakpoints,
    const double *slopes, Variable var) {
  assert(num_breakpoints >= 1);
  ++asl_->i.plterms_;
  plterm *term = Allocate<plterm>(
      sizeof(plterm) + 2 * num_breakpoints * sizeof(double));
  term->n = num_breakpoints + 1;
  double *data = term->bs;
  for (int i = 0; i < num_breakpoints; ++i) {
    data[2 * i] = slopes[i];
    data[2 * i + 1] = breakpoints[i];
  }
  data[2 * num_breakpoints] = slopes[num_breakpoints];
  ::expr *result = Allocate< ::expr>();
  result->op = r_ops_[OPPLTERM];
  result->L.p = term;
  result->R.e = var.expr_;
  return Expr::Create<PiecewiseLinearExpr>(result);
}

CallExpr ASLBuilder::MakeCall(Function f, ArrayRef<Expr> args) {
  int num_func_args = f.num_args();
  int num_args = SafeInt<int>(args.size()).value();
  if ((num_func_args >= 0 && num_args != num_func_args) ||
      (num_func_args < 0 && num_args < -(num_func_args + 1))) {
    throw Error("function {}: invalid number of arguments", f.name());
  }
  expr_f *result = Allocate<expr_f>(
        sizeof(expr_f) + SafeInt<int>(num_args - 1) * sizeof(::expr*));
  result->op = r_ops_[OPFUNCALL];
  result->fi = f.fi_;
  ::expr **arg_ptrs = result->args;
  int num_symbolic_args = 0, num_ifsyms = 0, num_constants = 0;
  for (int i = 0; i < num_args; ++i) {
    arg_ptrs[i] = args[i].expr_;
    if (Is<StringLiteral>(args[i]))
      ++num_symbolic_args;
    else if (args[i].kind() == expr::IFSYM)
      ++num_ifsyms;
    else if (Is<NumericConstant>(args[i]))
      ++num_constants;
  }
  num_symbolic_args += num_ifsyms;
  if (num_symbolic_args != 0 && (f.fi_->ftype & func::SYMBOLIC) == 0)
    throw Error("function {}: symbolic arguments not allowed", f.name());
  int num_numeric_args = num_args - num_symbolic_args;
  int kd = 0;
  double *ra = Allocate<double>(
      sizeof(arglist) +
      SafeInt<int>(num_constants + num_ifsyms) * sizeof(argpair) +
      SafeInt<int>(num_numeric_args + kd) * sizeof(double) +
      SafeInt<int>(num_symbolic_args) * sizeof(char*) +
      SafeInt<int>(num_args) * sizeof(int));
  double *b = ra + kd;
  arglist *al = result->al = reinterpret_cast<arglist*>(b + num_numeric_args);
  al->n = num_args;
  al->ra = ra;
  al->funcinfo = f.fi_->funcinfo;
  return Expr::Create<CallExpr>(reinterpret_cast< ::expr*>(result));
}

VarArgExpr ASLBuilder::MakeVarArg(expr::Kind kind, ArrayRef<NumericExpr> args) {
  CheckKind<VarArgExpr>(kind, "vararg");
  expr_va *result = Allocate<expr_va>();
  result->op = r_ops_[opcode(kind)];
  int num_args = SafeInt<int>(args.size()).value();
  de *d = result->L.d = Allocate<de>(
        SafeInt<int>(num_args) * sizeof(de) + sizeof(::expr*));
  for (int i = 0; i < num_args; ++i)
    d[i].e = args[i].expr_;
  d[num_args].e = 0;
  return Expr::Create<VarArgExpr>(reinterpret_cast< ::expr*>(result));
}

IteratedLogicalExpr ASLBuilder::MakeIteratedLogical(
    expr::Kind kind, ArrayRef<LogicalExpr> args) {
  CheckKind<IteratedLogicalExpr>(kind, "iterated logical");
  return MakeIterated<IteratedLogicalExpr>(kind, args);
}

StringLiteral ASLBuilder::MakeStringLiteral(fmt::StringRef value) {
  std::size_t size = value.size();
  expr_h *result = Allocate<expr_h>(
        SafeInt<int>(AddPadding(sizeof(expr_h)) + size));
  result->op = r_ops_[OPHOL];
  // Passing result->sym to std::copy causes an assertion failure in MSVC.
  char *dest = result->sym;
  const char *str = value.c_str();
  std::copy(str, str + size, dest);
  result->sym[size] = 0;
  return Expr::Create<StringLiteral>(reinterpret_cast< ::expr*>(result));
}
}
}  // namespace mp
