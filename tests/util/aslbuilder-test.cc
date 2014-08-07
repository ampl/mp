/*
 Tests of the ASL problem builder.

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

#include "gtest/gtest.h"

#include "solvers/util/aslbuilder.h"
#include "solvers/util/nl.h"
#include "tests/util.h"

#include <climits>

using ampl::NLHeader;
using ampl::Function;
using ampl::LogicalExpr;
using ampl::NumericExpr;
using ampl::MakeArrayRef;
using ampl::internal::ASLBuilder;
namespace func = ampl::func;

bool operator==(const cde &lhs, const cde &rhs) {
  return lhs.e == rhs.e && lhs.d == rhs.d && lhs.zaplen == rhs.zaplen;
}

bool operator==(const expr_v &lhs, const expr_v &rhs) {
  return lhs.op == rhs.op && lhs.a == rhs.a && lhs.v == rhs.v;
}

bool operator==(const cexp &lhs, const cexp &rhs) {
  return lhs.e == rhs.e && lhs.nlin == rhs.nlin && lhs.L == rhs.L &&
      lhs.funneled == rhs.funneled && lhs.cref == rhs.cref &&
      lhs.z.e == rhs.z.e && lhs.zlen == rhs.zlen && lhs.d == rhs.d &&
      lhs.vref == rhs.vref;
}

bool operator==(const cexp1 &lhs, const cexp1 &rhs) {
  return lhs.e == rhs.e && lhs.nlin == rhs.nlin && lhs.L == rhs.L;
}

namespace {

// Searches the ASL list for asl backwards from prev and forward from next
// returning true if it is found, false otherwise.
bool FindASL(ASLhead *prev, ASLhead *next, const ASL &asl) {
  while (prev && prev != &asl.p.h)
    prev = prev->prev;
  if (prev)
    return true;
  while (next && next != &asl.p.h)
    next = next->next;
  if (next)
    return true;
  return false;
}

std::size_t CountBlocks(void *start) {
  std::size_t num_blocks = 0;
  struct Mblock {
    struct Mblock *next;
    void *m[31];
  };
  for (Mblock *mb = reinterpret_cast<Mblock*>(start); mb; mb = mb->next)
    ++num_blocks;
  return num_blocks;
}

template <typename T>
T *GetPtr(const ASL &asl, T *Edaginfo::*ptr) { return asl.i.*ptr; }

template <typename T>
T *GetPtr(const ASL &asl, T *Edag1info::*ptr) {
  return reinterpret_cast<const ASL_fg&>(asl).I.*ptr;
}

// Checks if the array pointed to by ptr in the actual ASL object is the
// same as the one in the expected object. In particular, they should have
// the same offset from start and have the same elements.
template <typename StartInfo, typename Info, typename StartT, typename T>
void CheckArray(const ASL &expected, const ASL &actual,
    StartT *StartInfo::*start, T *Info::*ptr, std::size_t size,
    const char *str) {
  // Get start pointers.
  const char *expected_start =
      reinterpret_cast<const char*>(GetPtr(expected, start));
  const char *actual_start =
      reinterpret_cast<const char*>(GetPtr(actual, start));

  // Get array pointers and compare.
  const T *expected_ptr = GetPtr(expected, ptr);
  const T *actual_ptr = GetPtr(actual, ptr);
  EXPECT_EQ(reinterpret_cast<const char*>(expected_ptr) - expected_start,
      reinterpret_cast<const char*>(actual_ptr) - actual_start) << str;
  if (expected_ptr) {
    for (std::size_t i = 0; i < size; ++i)
      EXPECT_EQ(expected_ptr[i], actual_ptr[i]) << str << ' ' << i;
  }
}

#define CHECK_ARRAY(expected, actual, ptr, size) \
  CheckArray(expected, actual, &Edag1info::var_e_, ptr, size, #ptr)

template <typename T>
void ExpectArrayEqual(const T *expected,
    const T *actual, std::size_t size, const char *str) {
  if (!expected) {
    EXPECT_EQ(expected, actual);
    return;
  }
  for (std::size_t i = 0; i < size; ++i)
    EXPECT_EQ(expected[i], actual[i]) << str << ' ' << i;
}

// Compare two arrays for equality.
#define EXPECT_ARRAY_EQ(expected, actual, size) \
  ExpectArrayEqual(expected, actual, size, #expected)

// Compare two ASL objects for equality.
void CheckASL(const ASL &expected, const ASL &actual, bool complete = true) {
  // Compare Edagpars.
  EXPECT_TRUE(FindASL(expected.p.h.prev, expected.p.h.next, actual));
  EXPECT_TRUE(FindASL(actual.p.h.prev, actual.p.h.next, expected));
  EXPECT_EQ(expected.p.hffactor, actual.p.hffactor);
  EXPECT_EQ(expected.p.FUNNEL_MIN_, actual.p.FUNNEL_MIN_);
  EXPECT_EQ(expected.p.maxfwd_, actual.p.maxfwd_);
  EXPECT_EQ(expected.p.need_funcadd_, actual.p.need_funcadd_);
  EXPECT_EQ(expected.p.vrefGulp_, actual.p.vrefGulp_);
  EXPECT_EQ(expected.p.want_derivs_, actual.p.want_derivs_);
  EXPECT_EQ(expected.p.ihd_limit_, actual.p.ihd_limit_);
  EXPECT_EQ(expected.p.solve_code_, actual.p.solve_code_);
  EXPECT_EQ(expected.p.Objval, actual.p.Objval);
  EXPECT_EQ(expected.p.Objval_nomap, actual.p.Objval_nomap);
  EXPECT_EQ(expected.p.Objgrd, actual.p.Objgrd);
  EXPECT_EQ(expected.p.Objgrd_nomap, actual.p.Objgrd_nomap);
  EXPECT_EQ(expected.p.Conval, actual.p.Conval);
  EXPECT_EQ(expected.p.Jacval, actual.p.Jacval);
  EXPECT_EQ(expected.p.Conival, actual.p.Conival);
  EXPECT_EQ(expected.p.Conival_nomap, actual.p.Conival_nomap);
  EXPECT_EQ(expected.p.Congrd, actual.p.Congrd);
  EXPECT_EQ(expected.p.Congrd_nomap, actual.p.Congrd_nomap);
  EXPECT_EQ(expected.p.Hvcomp, actual.p.Hvcomp);
  EXPECT_EQ(expected.p.Hvcomp_nomap, actual.p.Hvcomp_nomap);
  EXPECT_EQ(expected.p.Hvinit, actual.p.Hvinit);
  EXPECT_EQ(expected.p.Hvinit_nomap, actual.p.Hvinit_nomap);
  EXPECT_EQ(expected.p.Hesset, actual.p.Hesset);
  EXPECT_EQ(expected.p.Lconval, actual.p.Lconval);
  EXPECT_EQ(expected.p.Xknown, actual.p.Xknown);
  EXPECT_EQ(expected.p.Duthes, actual.p.Duthes);
  EXPECT_EQ(expected.p.Duthes_nomap, actual.p.Duthes_nomap);
  EXPECT_EQ(expected.p.Fulhes, actual.p.Fulhes);
  EXPECT_EQ(expected.p.Fulhes_nomap, actual.p.Fulhes_nomap);
  EXPECT_EQ(expected.p.Sphes, actual.p.Sphes);
  EXPECT_EQ(expected.p.Sphes_nomap, actual.p.Sphes_nomap);
  EXPECT_EQ(expected.p.Sphset, actual.p.Sphset);
  EXPECT_EQ(expected.p.Sphset_nomap, actual.p.Sphset_nomap);

  // Compare Edaginfo.
  EXPECT_EQ(expected.i.ASLtype, actual.i.ASLtype);
  EXPECT_EQ(expected.i.amplflag_, actual.i.amplflag_);
  EXPECT_EQ(expected.i.need_nl_, actual.i.need_nl_);
  if (expected.i.ASLtype != ASL_read_f)
    CHECK_ARRAY(expected, actual, &Edaginfo::funcs_, expected.i.nfunc_);
  else
    EXPECT_EQ(expected.i.funcs_, actual.i.funcs_);
  EXPECT_EQ(expected.i.funcsfirst_, actual.i.funcsfirst_);
  EXPECT_EQ(expected.i.funcslast_, actual.i.funcslast_);
  EXPECT_EQ(expected.i.xscanf_, actual.i.xscanf_);
  for (int i = 0; i < NFHASH; ++i)
    EXPECT_EQ(expected.i.fhash_[i], actual.i.fhash_[i]);

  EXPECT_ARRAY_EQ(expected.i.adjoints_, actual.i.adjoints_, expected.i.amax_);
  EXPECT_EQ(expected.i.adjoints_nv1_ - expected.i.adjoints_,
      actual.i.adjoints_nv1_ - actual.i.adjoints_);
  EXPECT_ARRAY_EQ(expected.i.LUrhs_, actual.i.LUrhs_,
      2 * (expected.i.n_con_ + expected.i.nsufext[ASL_Sufkind_con]));
  EXPECT_EQ(expected.i.Urhsx_, actual.i.Urhsx_);
  EXPECT_EQ(expected.i.X0_, actual.i.X0_);

  // The number of variables and variable suffixes.
  std::size_t num_vars_suf = complete ?
      2 * (expected.i.n_var_ + expected.i.nsufext[ASL_Sufkind_var]) : 0;
  EXPECT_ARRAY_EQ(expected.i.LUv_, actual.i.LUv_, num_vars_suf);
  EXPECT_EQ(expected.i.Uvx_, actual.i.Uvx_);
  EXPECT_ARRAY_EQ(expected.i.Lastx_, actual.i.Lastx_, num_vars_suf);
  EXPECT_EQ(expected.i.pi0_, actual.i.pi0_);

  CHECK_ARRAY(expected, actual, &Edaginfo::objtype_, expected.i.n_obj_);
  EXPECT_EQ(expected.i.havex0_, actual.i.havex0_);
  EXPECT_EQ(expected.i.havepi0_, actual.i.havepi0_);
  EXPECT_EQ(expected.i.A_vals_, actual.i.A_vals_);
  EXPECT_EQ(expected.i.A_rownos_, actual.i.A_rownos_);
  EXPECT_EQ(expected.i.A_colstarts_, actual.i.A_colstarts_);

  EXPECT_EQ(expected.i.Cgrad_, actual.i.Cgrad_);
  CHECK_ARRAY(expected, actual, &Edaginfo::Ograd_, expected.i.n_obj_);
  EXPECT_EQ(expected.i.Cgrad0, actual.i.Cgrad0);

  EXPECT_EQ(expected.i.Fortran_, actual.i.Fortran_);
  EXPECT_EQ(expected.i.amax_, actual.i.amax_);

  EXPECT_EQ(expected.i.c_vars_, actual.i.c_vars_);
  EXPECT_EQ(expected.i.comb_, actual.i.comb_);
  EXPECT_EQ(expected.i.combc_, actual.i.combc_);
  EXPECT_EQ(expected.i.comc1_, actual.i.comc1_);
  EXPECT_EQ(expected.i.comc_, actual.i.comc_);
  EXPECT_EQ(expected.i.como1_, actual.i.como1_);
  EXPECT_EQ(expected.i.como_, actual.i.como_);

  EXPECT_EQ(expected.i.lnc_, actual.i.lnc_);
  EXPECT_EQ(expected.i.nbv_, actual.i.nbv_);
  EXPECT_EQ(expected.i.niv_, actual.i.niv_);
  EXPECT_EQ(expected.i.nlc_, actual.i.nlc_);
  EXPECT_EQ(expected.i.n_eqn_, actual.i.n_eqn_);
  EXPECT_EQ(expected.i.n_cc_, actual.i.n_cc_);
  EXPECT_EQ(expected.i.nlcc_, actual.i.nlcc_);
  EXPECT_EQ(expected.i.ndcc_, actual.i.ndcc_);

  EXPECT_EQ(expected.i.nzlb_, actual.i.nzlb_);
  EXPECT_EQ(expected.i.nlnc_, actual.i.nlnc_);
  EXPECT_EQ(expected.i.nlo_, actual.i.nlo_);
  EXPECT_EQ(expected.i.nlvb_, actual.i.nlvb_);
  EXPECT_EQ(expected.i.nlvc_, actual.i.nlvc_);
  EXPECT_EQ(expected.i.nlvo_, actual.i.nlvo_);
  EXPECT_EQ(expected.i.nlvbi_, actual.i.nlvbi_);
  EXPECT_EQ(expected.i.nlvci_, actual.i.nlvci_);
  EXPECT_EQ(expected.i.nlvoi_, actual.i.nlvoi_);
  EXPECT_EQ(expected.i.nwv_, actual.i.nwv_);
  EXPECT_EQ(expected.i.nzc_, actual.i.nzc_);
  EXPECT_EQ(expected.i.nzo_, actual.i.nzo_);
  EXPECT_EQ(expected.i.n_var_, actual.i.n_var_);
  EXPECT_EQ(expected.i.n_con_, actual.i.n_con_);
  EXPECT_EQ(expected.i.n_obj_, actual.i.n_obj_);
  EXPECT_EQ(expected.i.n_prob, actual.i.n_prob);
  EXPECT_EQ(expected.i.n_lcon_, actual.i.n_lcon_);
  EXPECT_EQ(expected.i.flags, actual.i.flags);
  EXPECT_EQ(expected.i.n_conjac_[0], actual.i.n_conjac_[0]);
  EXPECT_EQ(expected.i.n_conjac_[1], actual.i.n_conjac_[1]);

  EXPECT_EQ(expected.i.nclcon_, actual.i.nclcon_);
  EXPECT_EQ(expected.i.ncom0_, actual.i.ncom0_);
  EXPECT_EQ(expected.i.ncom1_, actual.i.ncom1_);
  EXPECT_EQ(expected.i.nderps_, actual.i.nderps_);
  EXPECT_EQ(expected.i.nfunc_, actual.i.nfunc_);
  EXPECT_EQ(expected.i.nzjac_, actual.i.nzjac_);
  EXPECT_EQ(expected.i.o_vars_, actual.i.o_vars_);
  EXPECT_EQ(expected.i.want_deriv_, actual.i.want_deriv_);
  EXPECT_EQ(expected.i.x0kind_, actual.i.x0kind_);
  EXPECT_EQ(expected.i.rflags, actual.i.rflags);
  EXPECT_EQ(expected.i.x0len_, actual.i.x0len_);

  EXPECT_STREQ(expected.i.filename_, actual.i.filename_);
  EXPECT_EQ(
      expected.i.stub_end_ - expected.i.filename_,
      actual.i.stub_end_ - actual.i.filename_);
  EXPECT_EQ(expected.i.archan_, actual.i.archan_);
  EXPECT_EQ(expected.i.awchan_, actual.i.awchan_);
  EXPECT_EQ(expected.i.binary_nl_, actual.i.binary_nl_);
  EXPECT_EQ(expected.i.return_nofile_, actual.i.return_nofile_);
  EXPECT_EQ(expected.i.plterms_, actual.i.plterms_);
  EXPECT_EQ(expected.i.maxrownamelen_, actual.i.maxrownamelen_);
  EXPECT_EQ(expected.i.maxcolnamelen_, actual.i.maxcolnamelen_);
  EXPECT_EQ(expected.i.co_index_, actual.i.co_index_);
  EXPECT_EQ(expected.i.cv_index_, actual.i.cv_index_);
  // Edaginfo::err_jmp_ is ignored.
  EXPECT_EQ(expected.i.err_jmp1_, actual.i.err_jmp1_);
  for (int i = 0; i < ampl::MAX_NL_OPTIONS + 1; ++i)
    EXPECT_EQ(expected.i.ampl_options_[i], actual.i.ampl_options_[i]);
  EXPECT_EQ(expected.i.obj_no_, actual.i.obj_no_);
  EXPECT_EQ(expected.i.nranges_, actual.i.nranges_);
  EXPECT_EQ(expected.i.want_xpi0_, actual.i.want_xpi0_);

  EXPECT_EQ(expected.i.c_cexp1st_, actual.i.c_cexp1st_);
  EXPECT_EQ(expected.i.o_cexp1st_, actual.i.o_cexp1st_);

  EXPECT_EQ(expected.i.cvar_, actual.i.cvar_);
  EXPECT_EQ(expected.i.ccind1, actual.i.ccind1);
  EXPECT_EQ(expected.i.ccind2, actual.i.ccind2);

  EXPECT_EQ(expected.i.size_expr_n_, actual.i.size_expr_n_);
  EXPECT_EQ(expected.i.ampl_vbtol_, actual.i.ampl_vbtol_);

  EXPECT_EQ(expected.i.zaC_, actual.i.zaC_);
  EXPECT_EQ(expected.i.zac_, actual.i.zac_);
  EXPECT_EQ(expected.i.zao_, actual.i.zao_);

  EXPECT_EQ(expected.i.skip_int_derivs_, actual.i.skip_int_derivs_);

  EXPECT_EQ(expected.i.nsuffixes, actual.i.nsuffixes);
  for (int i = 0; i < 4; ++i) {
    EXPECT_EQ(expected.i.nsufext[i], actual.i.nsufext[i]);
    EXPECT_EQ(expected.i.nsuff[i], actual.i.nsuff[i]);
    EXPECT_EQ(expected.i.suffixes[i], actual.i.suffixes[i]);
  }

  if (!expected.i.zerograds_) {
    EXPECT_EQ(expected.i.zerograds_, actual.i.zerograds_);
  } else {
    const char *expected_start =
        reinterpret_cast<const char*>(expected.i.zerograds_);
    const char *actual_start =
        reinterpret_cast<const char*>(actual.i.zerograds_);
    for (int i = 0; i < expected.i.n_obj_; ++i) {
      EXPECT_EQ(
          reinterpret_cast<const char*>(expected.i.zerograds_[i]) -
          expected_start,
          reinterpret_cast<const char*>(actual.i.zerograds_[i]) -
          actual_start);

      // Compare values pointed to by zerograds_[i].
      int nx = expected.i.nsufext[ASL_Sufkind_var];
      int nv = expected.i.nlvog;
      if (nv == 0) {
        nv = expected.i.n_var_;
        if (nv > expected.i.n_var0)
          nx -= nv - expected.i.n_var0;
      }
      int size = expected.i.n_obj_;
      for (ograd **ogp = expected.i.Ograd_, **ogpe = ogp + size;
          ogp < ogpe; ++ogp) {
        ograd *og = *ogp;
        int n = 0;
        while (og) {
          size += og->varno - n;
          n = og->varno + 1;
          if (n >= nv)
            break;
          og = og->next;
        }
        if (n < nv)
          size += nv - n;
      }
      size += expected.i.n_obj_ * nx;
      EXPECT_ARRAY_EQ(reinterpret_cast<const int*>(
          expected.i.zerograds_ + expected.i.n_obj_),
          reinterpret_cast<const int*>(actual.i.zerograds_ + expected.i.n_obj_),
          size);
    }
  }
  EXPECT_EQ(expected.i.congrd_mode, actual.i.congrd_mode);
  EXPECT_EQ(expected.i.x_known, actual.i.x_known);
  EXPECT_EQ(expected.i.xknown_ignore, actual.i.xknown_ignore);
  EXPECT_EQ(expected.i.zap_J, actual.i.zap_J);
  EXPECT_EQ(expected.i.nxval, actual.i.nxval);
  EXPECT_EQ(expected.i.nlvog, actual.i.nlvog);
  EXPECT_EQ(expected.i.ncxval, actual.i.ncxval);
  EXPECT_EQ(expected.i.noxval, actual.i.noxval);
  EXPECT_EQ(expected.i.sputinfo_, actual.i.sputinfo_);

  EXPECT_EQ(expected.i.Mblast - expected.i.Mbnext,
            actual.i.Mblast - actual.i.Mbnext);
  EXPECT_EQ(CountBlocks(expected.i.Mb), CountBlocks(actual.i.Mb));
  EXPECT_EQ(expected.i.memLast - expected.i.memNext,
      actual.i.memLast - actual.i.memNext);
  EXPECT_EQ(expected.i.ae, actual.i.ae);

  EXPECT_EQ(expected.i.connames, actual.i.connames);
  EXPECT_EQ(expected.i.lconnames, actual.i.lconnames);
  EXPECT_EQ(expected.i.objnames, actual.i.objnames);
  EXPECT_EQ(expected.i.varnames, actual.i.varnames);
  EXPECT_EQ(expected.i.vcochecked, actual.i.vcochecked);

  EXPECT_EQ(expected.i.uinfo, actual.i.uinfo);

  EXPECT_EQ(expected.i.iadjfcn, actual.i.iadjfcn);
  EXPECT_EQ(expected.i.dadjfcn, actual.i.dadjfcn);

  EXPECT_EQ(expected.i.cscale, actual.i.cscale);
  EXPECT_EQ(expected.i.vscale, actual.i.vscale);
  EXPECT_EQ(expected.i.lscale, actual.i.lscale);

  EXPECT_EQ(expected.i.arlast, actual.i.arlast);
  EXPECT_EQ(expected.i.arnext, actual.i.arnext);
  EXPECT_EQ(expected.i.arprev, actual.i.arprev);

  EXPECT_EQ(expected.i.csd, actual.i.csd);
  EXPECT_EQ(expected.i.rsd, actual.i.rsd);
  EXPECT_EQ(expected.i.n_var0, actual.i.n_var0);
  EXPECT_EQ(expected.i.n_con0, actual.i.n_con0);
  EXPECT_EQ(expected.i.n_var1, actual.i.n_var1);
  EXPECT_EQ(expected.i.n_con1, actual.i.n_con1);
  EXPECT_EQ(expected.i.vmap, actual.i.vmap);
  EXPECT_EQ(expected.i.cmap, actual.i.cmap);
  EXPECT_EQ(expected.i.vzap, actual.i.vzap);
  EXPECT_EQ(expected.i.czap, actual.i.czap);
  EXPECT_EQ(expected.i.vminv, actual.i.vminv);

  EXPECT_EQ(expected.i.Or, actual.i.Or);
  EXPECT_EQ(expected.i.orscratch, actual.i.orscratch);

  EXPECT_EQ(expected.i.mpa, actual.i.mpa);

  EXPECT_EQ(expected.i.Derrs, actual.i.Derrs);
  EXPECT_EQ(expected.i.Derrs0, actual.i.Derrs0);

  // Compare Edag1info.
  const ASL_fg &expected_fg = reinterpret_cast<const ASL_fg&>(expected);
  const ASL_fg &actual_fg = reinterpret_cast<const ASL_fg&>(actual);
  CHECK_ARRAY(expected, actual, &Edag1info::con_de_,
      expected.i.n_con_ + expected.i.nsufext[ASL_Sufkind_con]);
  CHECK_ARRAY(expected, actual, &Edag1info::lcon_de_, expected.i.n_lcon_);
  CHECK_ARRAY(expected, actual, &Edag1info::obj_de_, expected.i.n_obj_);
  CHECK_ARRAY(expected, actual, &Edag1info::var_e_,
      expected.i.n_var_ + expected.i.nsufext[ASL_Sufkind_var]);

  EXPECT_EQ(expected_fg.I.f_b_, actual_fg.I.f_b_);
  EXPECT_EQ(expected_fg.I.f_c_, actual_fg.I.f_c_);
  EXPECT_EQ(expected_fg.I.f_o_, actual_fg.I.f_o_);

  if (expected.i.ASLtype != ASL_read_f) {
    CHECK_ARRAY(expected, actual, &Edag1info::var_ex_, expected.i.ncom0_);
    CHECK_ARRAY(expected, actual, &Edag1info::var_ex1_,
        expected.i.comc1_ + expected.i.como1_);
    CHECK_ARRAY(expected, actual, &Edag1info::cexps_, expected.i.ncom0_);
    CHECK_ARRAY(expected, actual, &Edag1info::cexps1_, expected.i.ncom1_);
  } else {
    EXPECT_EQ(expected_fg.I.var_ex_, actual_fg.I.var_ex_);
    EXPECT_EQ(expected_fg.I.var_ex1_, actual_fg.I.var_ex1_);
    EXPECT_EQ(expected_fg.I.cexps_, actual_fg.I.cexps_);
    EXPECT_EQ(expected_fg.I.cexps1_, actual_fg.I.cexps1_);
  }

  EXPECT_EQ(expected_fg.I.r_ops_, actual_fg.I.r_ops_);
  EXPECT_EQ(expected_fg.I.c_class, actual_fg.I.c_class);
  EXPECT_EQ(expected_fg.I.o_class, actual_fg.I.o_class);
  EXPECT_EQ(expected_fg.I.v_class, actual_fg.I.v_class);
  EXPECT_EQ(expected_fg.I.c_class_max, actual_fg.I.c_class_max);
  EXPECT_EQ(expected_fg.I.o_class_max, actual_fg.I.o_class_max);
}

std::string HeaderToStr(const NLHeader &h) {
  fmt::Writer w;
  w << h;
  return w.str();
}

void CheckHeader(const NLHeader &h) {
  NLHeader actual_header = {};
  std::string nl = HeaderToStr(h);
  ampl::TextReader(nl, "(input)").ReadHeader(actual_header);

  EXPECT_EQ(h.format, actual_header.format);

  EXPECT_EQ(h.num_options, actual_header.num_options);
  for (int i = 0; i < ampl::MAX_NL_OPTIONS; ++i)
    EXPECT_EQ(h.options[i], actual_header.options[i]);
  EXPECT_EQ(h.ampl_vbtol, actual_header.ampl_vbtol);

  EXPECT_EQ(h.num_vars, actual_header.num_vars);
  EXPECT_EQ(h.num_algebraic_cons, actual_header.num_algebraic_cons);
  EXPECT_EQ(h.num_objs, actual_header.num_objs);
  EXPECT_EQ(h.num_ranges, actual_header.num_ranges);
  EXPECT_EQ(h.num_eqns, actual_header.num_eqns);
  EXPECT_EQ(h.num_logical_cons, actual_header.num_logical_cons);

  EXPECT_EQ(h.num_nl_cons, actual_header.num_nl_cons);
  EXPECT_EQ(h.num_nl_objs, actual_header.num_nl_objs);
  EXPECT_EQ(h.num_compl_conds, actual_header.num_compl_conds);
  EXPECT_EQ(h.num_nl_compl_conds, actual_header.num_nl_compl_conds);
  EXPECT_EQ(h.num_compl_dbl_ineqs, actual_header.num_compl_dbl_ineqs);
  EXPECT_EQ(h.num_compl_vars_with_nz_lb,
      actual_header.num_compl_vars_with_nz_lb);

  EXPECT_EQ(h.num_nl_net_cons, actual_header.num_nl_net_cons);
  EXPECT_EQ(h.num_linear_net_cons, actual_header.num_linear_net_cons);

  EXPECT_EQ(h.num_nl_vars_in_cons, actual_header.num_nl_vars_in_cons);
  EXPECT_EQ(h.num_nl_vars_in_objs, actual_header.num_nl_vars_in_objs);
  EXPECT_EQ(h.num_nl_vars_in_both, actual_header.num_nl_vars_in_both);

  EXPECT_EQ(h.num_linear_net_vars, actual_header.num_linear_net_vars);
  EXPECT_EQ(h.num_funcs, actual_header.num_funcs);
  EXPECT_EQ(h.flags, actual_header.flags);

  EXPECT_EQ(h.num_linear_binary_vars, actual_header.num_linear_binary_vars);
  EXPECT_EQ(h.num_linear_integer_vars, actual_header.num_linear_integer_vars);
  EXPECT_EQ(h.num_nl_integer_vars_in_both,
      actual_header.num_nl_integer_vars_in_both);
  EXPECT_EQ(h.num_nl_integer_vars_in_cons,
      actual_header.num_nl_integer_vars_in_cons);
  EXPECT_EQ(h.num_nl_integer_vars_in_objs,
      actual_header.num_nl_integer_vars_in_objs);

  EXPECT_EQ(h.num_con_nonzeros, actual_header.num_con_nonzeros);
  EXPECT_EQ(h.num_obj_nonzeros, actual_header.num_obj_nonzeros);

  EXPECT_EQ(h.max_con_name_len, actual_header.max_con_name_len);
  EXPECT_EQ(h.max_var_name_len, actual_header.max_var_name_len);

  EXPECT_EQ(h.num_common_exprs_in_both, actual_header.num_common_exprs_in_both);
  EXPECT_EQ(h.num_common_exprs_in_cons, actual_header.num_common_exprs_in_cons);
  EXPECT_EQ(h.num_common_exprs_in_objs, actual_header.num_common_exprs_in_objs);
  EXPECT_EQ(h.num_common_exprs_in_single_cons,
      actual_header.num_common_exprs_in_single_cons);
  EXPECT_EQ(h.num_common_exprs_in_single_objs,
      actual_header.num_common_exprs_in_single_objs);

  if (h.num_vars == 0)
    return;  // jac0dim fails if there are no vars

  WriteFile("test.nl", nl);
  char stub[] = "test.nl";
  ASL *asl = ASL_alloc(ASL_read_fg);
  jac0dim_ASL(asl, stub, static_cast<int>(strlen(stub)));
  std::remove(stub);

  EXPECT_EQ(asl->i.ampl_options_[0], actual_header.num_options);
  for (int i = 0; i < ampl::MAX_NL_OPTIONS; ++i)
    EXPECT_EQ(asl->i.ampl_options_[i + 1], actual_header.options[i]);
  EXPECT_EQ(asl->i.ampl_vbtol_, actual_header.ampl_vbtol);

  EXPECT_EQ(asl->i.n_var_, actual_header.num_vars);
  EXPECT_EQ(asl->i.n_con_, actual_header.num_algebraic_cons);
  EXPECT_EQ(asl->i.n_obj_, actual_header.num_objs);
  EXPECT_EQ(asl->i.nranges_, actual_header.num_ranges);
  EXPECT_EQ(asl->i.n_eqn_, actual_header.num_eqns);
  EXPECT_EQ(asl->i.n_lcon_, actual_header.num_logical_cons);

  EXPECT_EQ(asl->i.nlc_, actual_header.num_nl_cons);
  EXPECT_EQ(asl->i.nlo_, actual_header.num_nl_objs);
  EXPECT_EQ(asl->i.n_cc_, actual_header.num_compl_conds);
  EXPECT_EQ(asl->i.nlcc_, actual_header.num_nl_compl_conds);
  EXPECT_EQ(asl->i.ndcc_, actual_header.num_compl_dbl_ineqs);
  EXPECT_EQ(asl->i.nzlb_, actual_header.num_compl_vars_with_nz_lb);

  EXPECT_EQ(asl->i.nlnc_, actual_header.num_nl_net_cons);
  EXPECT_EQ(asl->i.lnc_, actual_header.num_linear_net_cons);

  EXPECT_EQ(asl->i.nlvc_, actual_header.num_nl_vars_in_cons);
  EXPECT_EQ(asl->i.nlvo_, actual_header.num_nl_vars_in_objs);
  EXPECT_EQ(asl->i.nlvb_, actual_header.num_nl_vars_in_both);

  EXPECT_EQ(asl->i.nwv_, actual_header.num_linear_net_vars);
  EXPECT_EQ(asl->i.nfunc_, actual_header.num_funcs);
  EXPECT_EQ(asl->i.flags, actual_header.flags);

  EXPECT_EQ(asl->i.nbv_, actual_header.num_linear_binary_vars);
  EXPECT_EQ(asl->i.niv_, actual_header.num_linear_integer_vars);
  EXPECT_EQ(asl->i.nlvbi_, actual_header.num_nl_integer_vars_in_both);
  EXPECT_EQ(asl->i.nlvci_, actual_header.num_nl_integer_vars_in_cons);
  EXPECT_EQ(asl->i.nlvoi_, actual_header.num_nl_integer_vars_in_objs);

  EXPECT_EQ(asl->i.nzc_, actual_header.num_con_nonzeros);
  EXPECT_EQ(asl->i.nzo_, actual_header.num_obj_nonzeros);

  EXPECT_EQ(asl->i.maxrownamelen_, actual_header.max_con_name_len);
  EXPECT_EQ(asl->i.maxcolnamelen_, actual_header.max_var_name_len);

  EXPECT_EQ(asl->i.comb_, actual_header.num_common_exprs_in_both);
  EXPECT_EQ(asl->i.comc_, actual_header.num_common_exprs_in_cons);
  EXPECT_EQ(asl->i.como_, actual_header.num_common_exprs_in_objs);
  EXPECT_EQ(asl->i.comc1_, actual_header.num_common_exprs_in_single_cons);
  EXPECT_EQ(asl->i.como1_, actual_header.num_common_exprs_in_single_objs);

  ASL_free(&asl);
}

TEST(NLTest, ReadFullHeader) {
  NLHeader header = {
    NLHeader::BINARY,
    9, {2, 3, 5, 7, 11, 13, 17, 19, 23}, 1.23,
    29, 47, 37, 41, 43, 31,
    53, 59, 67, 61, 71, 73,
    79, 83,
    89, 97, 101,
    103, 107, ampl::arith::GetKind(), 109,
    113, 127, 131, 137, 139,
    149, 151,
    157, 163,
    167, 173, 179, 181, 191
  };
  CheckHeader(header);
  NLHeader zero_header = {};
  CheckHeader(zero_header);
}

class ASLPtr {
 private:
  ASL *asl_;
  FMT_DISALLOW_COPY_AND_ASSIGN(ASLPtr);

 public:
  explicit ASLPtr(int type = ASL_read_fg) : asl_(ASL_alloc(type)) {}
  ~ASLPtr() { ASL_free(&asl_); }

  ASL *get() const { return asl_; }
  ASL &operator*() const { return *asl_; }
  ASL *operator->() const { return asl_; }
};

// Reads ASL header from a file with the specified header and body.
FILE *ReadHeader(ASL &asl, const NLHeader &h, const char *body) {
  fmt::Writer w;
  w << h << body;
  WriteFile("test.nl", w.str());
  char stub[] = "test";
  return jac0dim_ASL(&asl, stub, sizeof(stub) - 1);
}

// Check that ASLBuilder creates an ASL object compatible with the one
// created by jac0dim.
void CheckInitASL(const NLHeader &h) {
  ASLPtr expected, actual;
  fclose(ReadHeader(*expected, h, ""));
  ASLBuilder(actual.get()).InitASL("test", h);
  CheckASL(*expected, *actual);
}

TEST(ASLBuilderTest, Ctor) {
  ASLPtr asl;
  ASLBuilder builder(asl.get());
}

TEST(ASLBuilderTest, InitASLTrivial) {
  NLHeader header = {};
  header.num_vars = 1;  // jac0dim can't handle problems with 0 vars
  CheckInitASL(header);
}

TEST(ASLBuilderTest, InitASLFull) {
  NLHeader header = {
    NLHeader::BINARY,
    9, {2, 3, 5, 7, 11, 13, 17, 19, 23}, 1.23,
    29, 47, 37, 41, 43, 31,
    53, 59, 67, 61, 71, 73,
    79, 83,
    89, 97, 101,
    103, 107, ampl::arith::IEEE_BIG_ENDIAN, 109,
    113, 127, 131, 137, 139,
    149, 151,
    157, 163,
    167, 173, 179, 181, 191
  };
  CheckInitASL(header);
}

// Check that iadjfcn & dadjfcn are set properly when using different
// endianness.
TEST(ASLBuilderTest, ASLBuilderAdjFcn) {
  namespace arith = ampl::arith;
  arith::Kind arith_kind = arith::GetKind();
  if (arith::IsIEEE(arith_kind)) {
    NLHeader header = {NLHeader::BINARY};
    header.num_vars = 1;
    header.arith_kind = arith_kind == arith::IEEE_BIG_ENDIAN ?
          arith::IEEE_LITTLE_ENDIAN : arith::IEEE_BIG_ENDIAN;
    CheckInitASL(header);  // iadjfcn & dadjfcn are checked here.
  }
}

#define CHECK_THROW_ASL_ERROR(code, expected_error_code, expected_message) { \
  ampl::internal::ASLError error(0, ""); \
  try { \
    code; \
  } catch (const ampl::internal::ASLError &e) { \
    error = e; \
  } \
  EXPECT_EQ(expected_error_code, error.error_code()); \
  EXPECT_STREQ(expected_message, error.what()); \
}

TEST(ASLBuilderTest, ASLBuilderInvalidProblemDim) {
  NLHeader header = {};
  CHECK_THROW_ASL_ERROR(ASLBuilder().InitASL("test", header),
      ASL_readerr_corrupt, "invalid problem dimensions: M = 0, N = 0, NO = 0");
  header.num_vars = 1;
  ASLBuilder().InitASL("test", header);
  header.num_algebraic_cons = -1;
  CHECK_THROW_ASL_ERROR(ASLBuilder().InitASL("test", header),
      ASL_readerr_corrupt, "invalid problem dimensions: M = -1, N = 1, NO = 0");
  header.num_objs = -1;
  header.num_algebraic_cons = 0;
  CHECK_THROW_ASL_ERROR(ASLBuilder().InitASL("test", header),
      ASL_readerr_corrupt, "invalid problem dimensions: M = 0, N = 1, NO = -1");
}

// Check that x0len_ is set properly for different values of
// num_nl_vars_in_cons & num_nl_vars_in_objs.
TEST(ASLBuilderTest, ASLBuilderX0Len) {
  NLHeader header = {};
  header.num_vars = 1;
  header.num_nl_vars_in_cons = 5;
  header.num_nl_vars_in_objs = 10;
  CheckInitASL(header);
  std::swap(header.num_nl_vars_in_cons, header.num_nl_vars_in_objs);
  CheckInitASL(header);
}

int ReadASL(ASL &asl, const NLHeader &h, const char *body, int flags) {
  return fg_read_ASL(&asl, ReadHeader(asl, h, body), flags);
}

NLHeader MakeHeader(int num_vars = 1) {
  NLHeader header = {};
  header.num_vars = num_vars;
  header.num_objs = 1;
  return header;
}

TEST(ASLBuilderTest, ASLBuilderLinear) {
  NLHeader header = MakeHeader();
  ASLPtr actual(ASL_read_f);
  ASLBuilder builder(actual.get());
  builder.BeginBuild("test", header, 0);
  builder.EndBuild();
  ASLPtr expected(ASL_read_f);
  EXPECT_EQ(0,
      f_read_ASL(expected.get(), ReadHeader(*expected, header, ""), 0));
  CheckASL(*expected, *actual, false);
}

TEST(ASLBuilderTest, ASLBuilderTrivialProblem) {
  NLHeader header = MakeHeader();
  ASLPtr actual;
  ASLBuilder builder(actual.get());
  builder.BeginBuild("test", header, 0);
  builder.EndBuild();
  ASLPtr expected;
  EXPECT_EQ(0, ReadASL(*expected, header, "", 0));
  CheckASL(*expected, *actual, false);
}

TEST(ASLBuilderTest, ASLBuilderDisallowCLPByDefault) {
  NLHeader header = MakeHeader();
  header.num_logical_cons = 1;
  ASLPtr actual;
  ASLBuilder builder(actual.get());
  CHECK_THROW_ASL_ERROR(builder.BeginBuild("test", header, ASL_return_read_err),
      ASL_readerr_CLP, "cannot handle logical constraints");
  ASLPtr expected;
  EXPECT_EQ(ASL_readerr_CLP,
      ReadASL(*expected, header, "", ASL_return_read_err));
  CheckASL(*expected, *actual, false);
}

TEST(ASLBuilderTest, ASLBuilderAllowCLP) {
  NLHeader header = MakeHeader();
  header.num_logical_cons = 1;
  ASLPtr actual;
  ASLBuilder builder(actual.get());
  ampl::internal::ASLError error(0, "");
  builder.BeginBuild("test", header, ASL_return_read_err | ASL_allow_CLP);
  builder.EndBuild();
  ASLPtr expected;
  ReadASL(*expected, header, "", ASL_return_read_err | ASL_allow_CLP);
  CheckASL(*expected, *actual, false);
}

class TestASLBuilder : public ASLBuilder {
 public:
  enum {NUM_FUNCS = 2, NUM_OBJS = 3, NUM_CONS = 4, NUM_LOGICAL_CONS = 5};
  explicit TestASLBuilder(ASLPtr &asl, bool allow_missing_funcs = false)
  : ASLBuilder(asl.get()) {
    NLHeader header = MakeHeader();
    header.num_funcs = NUM_FUNCS;
    header.num_objs = NUM_OBJS;
    header.num_algebraic_cons = NUM_CONS;
    header.num_logical_cons = NUM_LOGICAL_CONS;
    int flags = ampl::internal::ASL_STANDARD_OPCODES | ASL_allow_CLP;
    if (allow_missing_funcs)
      flags |= ASL_allow_missing_funcs;
    BeginBuild("", header, flags);
  }
};

double TestFunc(arglist *) { return 0; }

template <std::size_t SIZE>
void CheckFunctionList(const ASLPtr &asl, const func_info *(&funcs)[SIZE]) {
  int index = 0;
  for (func_info *fi = asl->i.funcsfirst_; fi; fi = fi->fnext, ++index) {
    ASSERT_LT(index, SIZE);
    EXPECT_EQ(fi, funcs[index]);
  }
  EXPECT_EQ(SIZE, index);
}

TEST(ASLBuilderTest, SetObj) {
  ASLPtr asl;
  TestASLBuilder builder(asl);
#undef obj_de
  cde *obj_de = reinterpret_cast<ASL_fg*>(asl.get())->I.obj_de_;
  EXPECT_EQ(ampl::obj::MIN, asl->i.objtype_[1]);
  EXPECT_EQ(0, obj_de[1].e);
  builder.SetObj(1, ampl::obj::MAX, builder.MakeNumericConstant(42));
  EXPECT_EQ(ampl::obj::MAX, asl->i.objtype_[1]);
  EXPECT_EQ(reinterpret_cast<efunc*>(OPNUM), obj_de[1].e->op);
}

#ifndef NDEBUG
TEST(ASLBuilderTest, SetObjIndexOutOfRange) {
  ASLPtr asl;
  TestASLBuilder builder(asl);
  NumericExpr expr = builder.MakeNumericConstant(42);
  EXPECT_DEBUG_DEATH(builder.SetObj(-1, ampl::obj::MIN, expr), "Assertion");
  EXPECT_DEBUG_DEATH(builder.SetObj(TestASLBuilder::NUM_OBJS,
                                    ampl::obj::MIN, expr), "Assertion");
}
#endif

TEST(ASLBuilderTest, SetCon) {
  ASLPtr asl;
  TestASLBuilder builder(asl);
#undef con_de
  cde *con_de = reinterpret_cast<ASL_fg*>(asl.get())->I.con_de_;
  EXPECT_EQ(0, con_de[2].e);
  builder.SetCon(2, builder.MakeNumericConstant(42));
  EXPECT_EQ(reinterpret_cast<efunc*>(OPNUM), con_de[2].e->op);
}

#ifndef NDEBUG
TEST(ASLBuilderTest, SetConIndexOutOfRange) {
  ASLPtr asl;
  TestASLBuilder builder(asl);
  NumericExpr expr = builder.MakeNumericConstant(42);
  EXPECT_DEBUG_DEATH(builder.SetCon(-1, expr), "Assertion");
  EXPECT_DEBUG_DEATH(
        builder.SetCon(TestASLBuilder::NUM_CONS, expr), "Assertion");
}
#endif

TEST(ASLBuilderTest, SetLogicalCon) {
  ASLPtr asl;
  TestASLBuilder builder(asl);
#undef lcon_de
  cde *lcon_de = reinterpret_cast<ASL_fg*>(asl.get())->I.lcon_de_;
  EXPECT_EQ(0, lcon_de[2].e);
  builder.SetLogicalCon(2, builder.MakeLogicalConstant(true));
  EXPECT_EQ(reinterpret_cast<efunc*>(OPNUM), lcon_de[2].e->op);
}

#ifndef NDEBUG
TEST(ASLBuilderTest, SetLogicalConIndexOutOfRange) {
  ASLPtr asl;
  TestASLBuilder builder(asl);
  LogicalExpr expr = builder.MakeLogicalConstant(true);
  EXPECT_DEBUG_DEATH(builder.SetLogicalCon(-1, expr), "Assertion");
  EXPECT_DEBUG_DEATH(builder.SetLogicalCon(
                       TestASLBuilder::NUM_LOGICAL_CONS, expr), "Assertion");
}
#endif

TEST(ASLBuilderTest, AddFunction) {
  ASLPtr asl;
  TestASLBuilder builder(asl);
  char info = 0;
  EXPECT_EQ(0, asl->i.funcsfirst_);
  Function f = builder.AddFunction("foo", TestFunc, 11, func::SYMBOLIC, &info);
  EXPECT_EQ(11, f.num_args());
  EXPECT_STREQ("foo", f.name());
  builder.SetFunction(0, "foo", 11);
  func_info *fi = asl->i.funcs_[0];
  EXPECT_TRUE(TestFunc == fi->funcp);
  EXPECT_EQ(func::SYMBOLIC, fi->ftype);
  EXPECT_EQ(11, fi->nargs);
  EXPECT_EQ(&info, fi->funcinfo);
  const func_info *funcs[] = {fi};
  CheckFunctionList(asl, funcs);
  builder.AddFunction("bar", TestFunc, 0);
  builder.SetFunction(0, "bar", 0);
  const func_info *funcs2[] = {fi, asl->i.funcs_[0]};
  CheckFunctionList(asl, funcs2);
  EXPECT_NE(Function(), builder.AddFunction("f1", TestFunc, 0, func::NUMERIC));
}

TEST(ASLBuilderTest, SetFunction) {
  ASLPtr asl;
  TestASLBuilder builder(asl);
  builder.AddFunction("foo", TestFunc, 3, func::SYMBOLIC);
  EXPECT_EQ(0, asl->i.funcs_[1]);
  Function f = builder.SetFunction(1, "foo", 3, func::SYMBOLIC);
  EXPECT_STREQ("foo", f.name());
  EXPECT_EQ(3, f.num_args());
  EXPECT_STREQ("foo", asl->i.funcs_[1]->name);
  EXPECT_EQ(3, asl->i.funcs_[1]->nargs);
  EXPECT_EQ(1, asl->i.funcs_[1]->ftype);
}

TEST(ASLBuilderTest, SetFunctionSameIndex) {
  ASLPtr asl;
  TestASLBuilder builder(asl);
  builder.AddFunction("f", TestFunc, 0);
  builder.AddFunction("g", TestFunc, 0);
  EXPECT_EQ(0, asl->i.funcs_[0]);
  builder.SetFunction(0, "f", 0);
  EXPECT_STREQ("f", asl->i.funcs_[0]->name);
  builder.SetFunction(0, "g", 0);
  EXPECT_STREQ("g", asl->i.funcs_[0]->name);
}

#ifndef NDEBUG
TEST(ASLBuilderTest, SetFunctionIndexOutOfRange) {
  ASLPtr asl;
  TestASLBuilder builder(asl);
  builder.AddFunction("foo", TestFunc, 3);
  EXPECT_DEBUG_DEATH(builder.SetFunction(-1, "f", 0), "Assertion");
  EXPECT_DEBUG_DEATH(
        builder.SetFunction(TestASLBuilder::NUM_FUNCS, "foo", 3), "Assertion");
}
#endif

TEST(ASLBuilderTest, SetFunctionMatchNumArgs) {
  ASLPtr asl;
  TestASLBuilder builder(asl);
  builder.AddFunction("f", TestFunc, 3);
  CHECK_THROW_ASL_ERROR(builder.SetFunction(0, "f", 2),
    ASL_readerr_argerr, "function f: disagreement of nargs: 3 and 2");
  builder.AddFunction("g", TestFunc, -3);
  builder.SetFunction(0, "g", 0);
  builder.SetFunction(0, "f", -4);
  CHECK_THROW_ASL_ERROR(builder.SetFunction(0, "f", -5),
    ASL_readerr_argerr, "function f: disagreement of nargs: 3 and -5");
}

TEST(ASLBuilderTest, SetMissingFunction) {
  ASLPtr asl;
  TestASLBuilder builder(asl);
  CHECK_THROW_ASL_ERROR(builder.SetFunction(0, "f", 0),
    ASL_readerr_unavail, "function f not available");
  ASLPtr asl2;
  TestASLBuilder builder2(asl2, true);
  builder2.SetFunction(0, "f", 0);
  EXPECT_STREQ("f", asl2->i.funcs_[0]->name);
  arglist al = {};
  asl2->i.funcs_[0]->funcp(&al);
  EXPECT_STREQ("attempt to call unavailable function", al.Errmsg);
}

// The ASLBuilder::Make* methods are tested in expr-test because they serve
// the purpose of expression constructors and testing them separately from
// expression classes doesn't make much sense.

TEST(ASLBuilderTest, SizeOverflow) {
  using safeint::OverflowError;
  ASLBuilder builder;
  NLHeader header = MakeHeader();
  header.num_funcs = 1;
  builder.BeginBuild("", header,
    ampl::internal::ASL_STANDARD_OPCODES | ASL_allow_missing_funcs);
  NumericExpr args[] = {builder.MakeNumericConstant(0)};

  std::size_t num_args = (INT_MAX - sizeof(expr_f)) / sizeof(expr*) + 2;
  std::size_t max_size = INT_MAX + 1u;
  Function f = builder.AddFunction("f", TestFunc, -1);
  EXPECT_THROW(builder.MakeCall(
                 f, MakeArrayRef<ampl::Expr>(args, num_args)), OverflowError);
  EXPECT_THROW(builder.MakeCall(
                 f, MakeArrayRef<ampl::Expr>(args, max_size)), OverflowError);

  num_args = (INT_MAX - sizeof(expr*)) / sizeof(de) + 1;
  EXPECT_THROW(builder.MakeVarArg(
                 MINLIST, MakeArrayRef(args, num_args)), OverflowError);
  EXPECT_THROW(builder.MakeVarArg(
                 MINLIST, MakeArrayRef(args, max_size)), OverflowError);

  num_args = (INT_MAX - sizeof(expr) + sizeof(double)) / sizeof(expr*);
  EXPECT_THROW(builder.MakeSum(MakeArrayRef(args, num_args)), OverflowError);
  EXPECT_THROW(builder.MakeSum(MakeArrayRef(args, max_size)), OverflowError);

  LogicalExpr largs[] = {builder.MakeLogicalConstant(false)};
  EXPECT_THROW(builder.MakeCount(MakeArrayRef(largs, num_args)), OverflowError);
  EXPECT_THROW(builder.MakeCount(MakeArrayRef(largs, max_size)), OverflowError);

  EXPECT_THROW(builder.MakeNumberOf(
                 MakeArrayRef(args, num_args)), OverflowError);
  EXPECT_THROW(builder.MakeNumberOf(
                 MakeArrayRef(args, max_size)), OverflowError);

  EXPECT_THROW(builder.MakeIteratedLogical(
                 ORLIST, MakeArrayRef(largs, num_args)), OverflowError);
  EXPECT_THROW(builder.MakeIteratedLogical(
                 ORLIST, MakeArrayRef(largs, max_size)), OverflowError);

  EXPECT_THROW(builder.MakeAllDiff(
                 MakeArrayRef(args, num_args)), OverflowError);
  EXPECT_THROW(builder.MakeAllDiff(
                 MakeArrayRef(args, max_size)), OverflowError);

  std::size_t size = INT_MAX - sizeof(expr_h) + 1;
  EXPECT_THROW(builder.MakeStringLiteral(
                 fmt::StringRef("", size)), OverflowError);
  EXPECT_THROW(builder.MakeStringLiteral(
                 fmt::StringRef("", max_size)), OverflowError);
}

// Test that ASLBuilder can act as a handler for NLReader.
TEST(ASLBuilderTest, NLHandler) {
  ASLBuilder builder;
  ampl::ReadNLString(HeaderToStr(MakeHeader()), builder);
}

// TODO: test SetVarBounds, SetConBounds, AddSuffix
}
