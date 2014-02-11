/*
 Tests of the the AMPL expression factory.

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

#include "solvers/util/expr-factory.h"
#include "solvers/util/nl.h"
#include "tests/util.h"

using ampl::ExprFactory;
using ampl::NLHeader;

bool operator==(const cde &lhs, const cde &rhs) {
  return lhs.e == rhs.e && lhs.d == rhs.d && lhs.zaplen == rhs.zaplen;
}

bool operator==(const expr_v &lhs, const expr_v &rhs) {
  return lhs.op == rhs.op && lhs.a == rhs.a && lhs.v == rhs.v;
}

namespace {

TEST(ExprFactoryTest, Ctor) {
  NLHeader h = {};
  h.num_objs = 1;
  ExprFactory ef(h, "");
}

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

void GetInfo(const ASL &asl, const Edaginfo *&info) { info = &asl.i; }
void GetInfo(const ASL &asl, const Edag1info *&info) {
  info = &reinterpret_cast<const ASL_fg&>(asl).I;
}

template <typename Info, typename T>
void CheckArray(const ASL &expected, const ASL &actual,
    T *Info::*ptr, std::size_t size, const char *str) {
  const char *expected_start = reinterpret_cast<const char*>(
      reinterpret_cast<const ASL_fg&>(expected).I.var_e_);
  const char *actual_start = reinterpret_cast<const char*>(
      reinterpret_cast<const ASL_fg&>(actual).I.var_e_);
  const Info *expected_info = 0, *actual_info = 0;
  GetInfo(expected, expected_info);
  GetInfo(actual, actual_info);
  const char *expected_ptr = reinterpret_cast<const char*>(expected_info->*ptr);
  const char *actual_ptr = reinterpret_cast<const char*>(actual_info->*ptr);
  EXPECT_EQ(expected_ptr - expected_start, actual_ptr - actual_start) << str;
  if (expected_ptr) {
    for (std::size_t i = 0; i < size; ++i) {
      EXPECT_EQ((expected_info->*ptr)[i], (actual_info->*ptr)[i])
          << str << ' ' << i;
    }
  }
}

#define CHECK_ARRAY(expected, actual, ptr, size) \
  CheckArray(expected, actual, ptr, size, #ptr)

template <typename T>
void ExpectArrayEqual(const T *expected,
    const T *actual, std::size_t size, const char *str) {
  if (!expected) {
    EXPECT_EQ(expected, actual);
    return;
  }
  for (std::size_t i = 0; i < size; ++i)
    EXPECT_EQ(expected[i], actual[i]) << str;
}

// Compare two arrays for equality.
#define EXPECT_ARRAY_EQ(expected, actual, size) \
  ExpectArrayEqual(expected, actual, size, #expected)

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
  CHECK_ARRAY(expected, actual, &Edaginfo::funcs_, expected.i.nfunc_);
  EXPECT_EQ(expected.i.funcsfirst_, actual.i.funcsfirst_);
  EXPECT_EQ(expected.i.funcslast_, actual.i.funcslast_);
  EXPECT_EQ(expected.i.xscanf_, actual.i.xscanf_);
  for (int i = 0; i < NFHASH; ++i)
    EXPECT_EQ(expected.i.fhash_[i], actual.i.fhash_[i]);

  EXPECT_ARRAY_EQ(expected.i.adjoints_, actual.i.adjoints_, expected.i.amax_);
  EXPECT_EQ(expected.i.adjoints_nv1_ - expected.i.adjoints_,
      actual.i.adjoints_nv1_ - actual.i.adjoints_);
  EXPECT_ARRAY_EQ(expected.i.LUrhs_, actual.i.LUrhs_, 2 * sizeof(double) *
      (expected.i.n_con_ + expected.i.nsufext[ASL_Sufkind_con]));
  EXPECT_EQ(expected.i.Urhsx_, actual.i.Urhsx_);
  EXPECT_EQ(expected.i.X0_, actual.i.X0_);
  std::size_t luv_size = complete ? 2 * sizeof(double) *
      (expected.i.n_var_ + expected.i.nsufext[ASL_Sufkind_var]) : 0;
  EXPECT_ARRAY_EQ(expected.i.LUv_, actual.i.LUv_, luv_size);
  EXPECT_EQ(expected.i.Uvx_, actual.i.Uvx_);
  EXPECT_EQ(expected.i.Lastx_, actual.i.Lastx_);
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
  EXPECT_EQ(expected.i.err_jmp_, actual.i.err_jmp_);
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

  EXPECT_EQ(expected.i.zerograds_, actual.i.zerograds_);
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
  EXPECT_EQ(expected.i.memNext, actual.i.memNext);
  EXPECT_EQ(expected.i.memLast, actual.i.memLast);
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
      expected.i.n_var_ + expected.i.nsufext[ASL_Sufkind_var] +
      expected.i.comb_ + expected.i.comc_ + expected.i.como_ +
      expected.i.comc1_ + expected.i.como1_);

  EXPECT_EQ(expected_fg.I.f_b_, actual_fg.I.f_b_);
  EXPECT_EQ(expected_fg.I.f_c_, actual_fg.I.f_c_);
  EXPECT_EQ(expected_fg.I.f_o_, actual_fg.I.f_o_);
  EXPECT_EQ(expected_fg.I.var_ex_, actual_fg.I.var_ex_);
  EXPECT_EQ(expected_fg.I.var_ex1_, actual_fg.I.var_ex1_);
  EXPECT_EQ(expected_fg.I.cexps_, actual_fg.I.cexps_);
  EXPECT_EQ(expected_fg.I.cexps1_, actual_fg.I.cexps1_);
  EXPECT_EQ(expected_fg.I.r_ops_, actual_fg.I.r_ops_);
  EXPECT_EQ(expected_fg.I.c_class, actual_fg.I.c_class);
  EXPECT_EQ(expected_fg.I.o_class, actual_fg.I.o_class);
  EXPECT_EQ(expected_fg.I.v_class, actual_fg.I.v_class);
  EXPECT_EQ(expected_fg.I.c_class_max, actual_fg.I.c_class_max);
  EXPECT_EQ(expected_fg.I.o_class_max, actual_fg.I.o_class_max);
}

class ASLPtr : ampl::Noncopyable {
 private:
  ASL *asl_;

 public:
  ASLPtr() : asl_(ASL_alloc(ASL_read_fg)) {}
  ~ASLPtr() { ASL_free(&asl_); }

  ASL *get() const { return asl_; }
  ASL &operator*() const { return *asl_; }
};

// Reads ASL header from a file with the specified header and body.
FILE *ReadHeader(ASL &asl, const NLHeader &h, const char *body) {
  fmt::Writer w;
  w << h << body;
  WriteFile("test.nl", w.str());
  char stub[] = "test";
  return jac0dim_ASL(&asl, stub, sizeof(stub) - 1);
}

// Check that InitASL creates an ASL object compatible with the one
// created by jac0dim.
void CheckInitASL(const NLHeader &h) {
  ASLPtr expected, actual;
  fclose(ReadHeader(*expected, h, ""));
  ampl::internal::InitASL(*actual, "test", h);
  CheckASL(*expected, *actual);
}

TEST(ExprFactoryTest, InitASLTrivial) {
  NLHeader header = {};
  header.num_vars = 1;  // jac0dim can't handle problems with 0 vars
  CheckInitASL(header);
}

TEST(ExprFactoryTest, InitASLFull) {
  NLHeader header = {
    NLHeader::BINARY,
    9, {2, 3, 5, 7, 11, 13, 17, 19, 23}, 1.23,
    29, 47, 37, 41, 43, 31,
    53, 59, 67, 61, 71, 73,
    79, 83,
    89, 97, 101,
    103, 107, 109,
    113, 127, 131, 137, 139,
    149, 151,
    157, 163,
    167, 173, 179, 181, 191
  };
  CheckInitASL(header);
}

TEST(ExprFactoryTest, InitASLAdjFcn) {
  NLHeader header = {NLHeader::BINARY_SWAPPED};
  header.num_vars = 1;
  CheckInitASL(header);
}

TEST(ExprFactoryTest, ASLBuilder) {
  NLHeader header = {};
  header.num_vars = header.num_objs = 1;
  ASLPtr actual;
  ampl::internal::InitASL(*actual, "test", header);
  ampl::internal::ASLBuilder builder(*actual);
  builder.BeginBuild(0);
  builder.EndBuild();
  ASLPtr expected;
  fg_read_ASL(expected.get(), ReadHeader(*expected, header, ""), 0);
  CheckASL(*expected, *actual, false);
}

TEST(ExprFactoryTest, CreateNumericConstant) {
  NLHeader header = {};
  header.num_vars = header.num_objs = 1;
  ExprFactory ef(header, "");
  EXPECT_EQ(42.0, ef.CreateNumericConstant(42).value());
}

TEST(ExprFactoryTest, CreateVariable) {
  NLHeader header = {};
  header.num_vars = 10;
  header.num_objs = 1;
  ExprFactory ef(header, "");
  EXPECT_EQ(0, ef.CreateVariable(0).index());
  EXPECT_EQ(9, ef.CreateVariable(9).index());
  EXPECT_DEBUG_DEATH(ef.CreateVariable(-1);, "Assertion");  // NOLINT(*)
  EXPECT_DEBUG_DEATH(ef.CreateVariable(10);, "Assertion");  // NOLINT(*)
}

// TODO: check if the asl produced by ExprFactory is binary compatible with
//       the one produced by reading an .nl file with fg_read.

// TODO: more tests
}
