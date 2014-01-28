/*
 An .nl reader tests.

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

#include "gtest/gtest.h"
#include "solvers/util/nlreader.h"
#include "solvers/util/problem.h"
#include "tests/util.h"

using ampl::NLHeader;
using ampl::NLReader;

namespace {

struct TestNLHandler : ampl::NLHandler {
  NLHeader header;
  void HandleHeader(const NLHeader &h) { header = h; }
};

std::string ReplaceLine(std::string s, int line_index, const char *new_line) {
  std::string::size_type start = 0;
  while (line_index-- > 0) {
    start = s.find('\n', start);
    if (start == std::string::npos)
      throw ampl::Error("invalid line index");
    ++start;
  }
  std::string::size_type end = s.find('\n', start);
  if (end == std::string::npos)
    end = s.size();
  s.replace(start, end - start, new_line);
  return s;
}

TEST(NLReaderTest, ReplaceLine) {
  EXPECT_EQ("de", ReplaceLine("", 0, "de"));
  EXPECT_EQ("de", ReplaceLine("abc", 0, "de"));
  EXPECT_THROW(ReplaceLine("abc", 1, "de"), ampl::Error);
  EXPECT_EQ("de\n", ReplaceLine("abc\n", 0, "de"));
  EXPECT_EQ("abc\nde", ReplaceLine("abc\n", 1, "de"));
  EXPECT_EQ("gh\ndef", ReplaceLine("abc\ndef", 0, "gh"));
  EXPECT_EQ("abc\ngh", ReplaceLine("abc\ndef", 1, "gh"));
  EXPECT_THROW(ReplaceLine("abc\ndef", 2, "gh"), ampl::Error);
}

// Formats header as a string.
std::string FormatHeader(const NLHeader &h) {
  fmt::Writer w;
  w << 'g' << h.num_options;
  for (int i = 0; i < ampl::MAX_NL_OPTIONS; ++i)
    w << ' ' << h.options[i];
  w << ' ' << h.ampl_vbtol << '\n';
  w.Format(" {} {} {} {} {} {}\n")
      << h.num_vars << (h.num_cons - h.num_logical_cons) << h.num_objs
      << h.num_ranges << h.num_eqns << h.num_logical_cons;
  w.Format(" {} {} {} {} {} {}\n")
      << h.num_nl_cons << h.num_nl_objs
      << (h.num_compl_conds - h.num_nl_compl_conds)
      << h.num_nl_compl_conds << h.num_compl_dbl_ineqs
      << h.num_compl_vars_with_nz_lb;
  w.Format(" {} {}\n")
      << h.num_nl_net_cons << h.num_linear_net_cons;
  w.Format(" {} {} {}\n")
      << h.num_nl_vars_in_cons << h.num_nl_vars_in_objs
      << h.num_nl_vars_in_both;
  w.Format(" {} {} 0 {}\n")
      << h.num_linear_net_vars << h.num_funcs << h.flags;
  w.Format(" {} {} {} {} {}\n")
      << h.num_linear_binary_vars << h.num_linear_integer_vars
      << h.num_nl_integer_vars_in_both << h.num_nl_integer_vars_in_cons
      << h.num_nl_integer_vars_in_objs;
  w.Format(" {} {}\n")
      << h.num_con_nonzeros << h.num_obj_nonzeros;
  w.Format(" {} {}\n")
      << h.max_con_name_len << h.max_var_name_len;
  w.Format(" {} {} {} {} {}\n")
      << h.num_common_exprs_in_both << h.num_common_exprs_in_cons
      << h.num_common_exprs_in_objs << h.num_common_exprs_in_cons1
      << h.num_common_exprs_in_objs1;
  return w.str();
}

// Reads a zero header with one modified line.
NLHeader ReadHeader(int line_index, const char *line) {
  TestNLHandler handler;
  NLReader reader(&handler);
  reader.ReadString(ReplaceLine(FormatHeader(NLHeader()), line_index, line));
  return handler.header;
}

TEST(NLReaderTest, NoNewlineAtEOF) {
  NLReader().ReadString("g\n"
    " 1 1 0\n"
    " 0 0\n"
    " 0 0\n"
    " 0 0 0\n"
    " 0 0 0 1\n"
    " 0 0 0 0 0\n"
    " 0 0\n"
    " 0 0\n"
    " 0 0 0 0 0\n"
    "k0\0h");
}

TEST(NLReaderTest, InvalidFormat) {
  EXPECT_THROW_MSG(ReadHeader(0, "x"),
      ampl::ParseError, "(input):1:1: invalid format 'x'");
}

TEST(NLReaderTest, InvalidNumOptions) {
  EXPECT_EQ(0, ReadHeader(0, "ga").num_options);
  EXPECT_EQ(0, ReadHeader(0, "g-1").num_options);
  EXPECT_THROW_MSG(ReadHeader(0, "g10"),
      ampl::ParseError, "(input):1:2: too many options");
  EXPECT_THROW_MSG(ReadHeader(0,
      c_str(fmt::Format("g{}") << static_cast<unsigned>(INT_MAX) + 1)),
      ampl::ParseError, "(input):1:2: number is too big");
}

void CheckReadOptions(size_t num_options,
    size_t num_options_to_write, const int *options) {
  fmt::Writer w;
  w << 'g' << num_options;
  for (size_t i = 0; i < num_options_to_write; ++i)
    w << ' ' << options[i];
  NLHeader header = ReadHeader(0, w.c_str());
  ASSERT_EQ(num_options, header.num_options);
  size_t min_num_options = std::min(num_options, num_options_to_write);
  for (size_t i = 0; i < min_num_options; ++i)
    EXPECT_EQ(options[i], header.options[i]);
  for (size_t i = min_num_options; i < num_options_to_write; ++i)
    EXPECT_EQ(0, header.options[i]);
}

TEST(NLReaderTest, ReadOptions) {
  const int options[ampl::MAX_NL_OPTIONS + 1] = {
      3, 5, 7, 11, 13, 17, 19, 23, 29, 31
  };
  for (size_t i = 0; i < ampl::MAX_NL_OPTIONS; ++i) {
    for (size_t j = 0; j < ampl::MAX_NL_OPTIONS + 1; ++j)
      CheckReadOptions(i, j, options);
  }
  EXPECT_EQ(0, ReadHeader(0, "g").num_options);
}

TEST(NLReaderTest, ReadAMPLVBTol) {
  EXPECT_EQ(4.2, ReadHeader(0, "g2 0 3 4.2").ampl_vbtol);
  EXPECT_EQ(0, ReadHeader(0, "g2 0 0 4.2").ampl_vbtol);
  EXPECT_EQ(0, ReadHeader(0, "g2 0 3").ampl_vbtol);
}

void CheckHeader(const NLHeader &h) {
  std::string nl = FormatHeader(h);
  TestNLHandler handler;
  NLReader reader(&handler);
  reader.ReadString(nl);
  NLHeader actual_header = handler.header;

  EXPECT_EQ(h.num_options, actual_header.num_options);
  for (int i = 0; i < ampl::MAX_NL_OPTIONS; ++i)
    EXPECT_EQ(h.options[i], actual_header.options[i]);
  EXPECT_EQ(h.ampl_vbtol, actual_header.ampl_vbtol);

  EXPECT_EQ(h.num_vars, actual_header.num_vars);
  EXPECT_EQ(h.num_cons, actual_header.num_cons);
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
  // TODO: test arith_kind
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
  EXPECT_EQ(h.num_common_exprs_in_cons1,
      actual_header.num_common_exprs_in_cons1);
  EXPECT_EQ(h.num_common_exprs_in_objs1,
      actual_header.num_common_exprs_in_objs1);

  if (h.num_vars == 0)
    return;  // jac0dim fails if there are no vars

  WriteFile("test.nl", nl);
  char stub[] = "test.nl";
  ASL *asl = ASL_alloc(ASL_read_fg);
  jac0dim_ASL(asl, stub, strlen(stub));
  std::remove(stub);

  EXPECT_EQ(asl->i.ampl_options_[0], actual_header.num_options);
  for (int i = 0; i < ampl::MAX_NL_OPTIONS; ++i)
    EXPECT_EQ(asl->i.ampl_options_[i + 1], actual_header.options[i]);
  EXPECT_EQ(asl->i.ampl_vbtol_, actual_header.ampl_vbtol);

  EXPECT_EQ(asl->i.n_var_, actual_header.num_vars);
  EXPECT_EQ(asl->i.n_con_,
      actual_header.num_cons - actual_header.num_logical_cons);
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
  // TODO: test arith_kind
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
  EXPECT_EQ(asl->i.comc1_, actual_header.num_common_exprs_in_cons1);
  EXPECT_EQ(asl->i.como1_, actual_header.num_common_exprs_in_objs1);

  ASL_free(&asl);
}

TEST(NLReaderTest, ReadFullHeader) {
  NLHeader header = {
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
  CheckHeader(header);
  NLHeader zero_header = {};
  CheckHeader(zero_header);
}

TEST(NLReaderTest, IncompleteHeader) {
  EXPECT_THROW_MSG(
      ReadHeader(1, " 1 0"),
      ampl::ParseError, "(input):2:5: expected nonnegative integer");
  EXPECT_THROW_MSG(
      ReadHeader(2, " 0"),
      ampl::ParseError, "(input):3:3: expected nonnegative integer");
  EXPECT_THROW_MSG(
      ReadHeader(3, " 0"),
      ampl::ParseError, "(input):4:3: expected nonnegative integer");
  EXPECT_THROW_MSG(
      ReadHeader(4, " 0"),
      ampl::ParseError, "(input):5:3: expected nonnegative integer");
  // TODO: test required and optional elements in all lines
  //       also check the defaults for optional elements
}

TEST(NLReaderTest, NumComplDblIneq) {
  EXPECT_EQ(42, ReadHeader(2, " 0 0 0 0 42").num_compl_dbl_ineqs);
  EXPECT_EQ(-1, ReadHeader(2, " 0 0 70 0 42").num_compl_dbl_ineqs);
}

// TODO: more tests
}
