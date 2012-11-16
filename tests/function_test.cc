/*
 Tests of the AMPL function testing infrastructure.

 Copyright (C) 2012 AMPL Optimization LLC

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization LLC disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#include <limits>
#include <sstream>
#include <cmath>
#include <cstring>

#include "gtest/gtest.h"
#include "tests/function.h"
#include "tests/config.h"
#include "solvers/asl.h"

using std::ptr_fun;
using std::sqrt;
using std::vector;

using fun::BitSet;
using fun::DerivativeBinder;
using fun::Differentiator;
using fun::Function;
using fun::FunctionInfo;
using fun::FunctionPointer;
using fun::FunctionWithTypes;
using fun::GetType;
using fun::Handler;
using fun::Library;
using fun::MakeArgs;
using fun::MakeVariant;
using fun::Tuple;
using fun::Variant;
using fun::Table;

namespace {

#undef VOID

TEST(FunctionTest, Type) {
  EXPECT_EQ(fun::VOID, GetType<void>::VALUE);
  EXPECT_EQ(fun::INT, GetType<int>::VALUE);
  EXPECT_EQ(fun::UINT, GetType<unsigned>::VALUE);
  EXPECT_EQ(fun::DOUBLE, GetType<double>::VALUE);
  EXPECT_EQ(fun::STRING, GetType<const char*>::VALUE);
  EXPECT_EQ(fun::POINTER, GetType<int*>::VALUE);
}

TEST(FunctionTest, VariantCtor) {
  Variant v;
  EXPECT_EQ(fun::VOID, v.type());
  EXPECT_THROW(v.number(), std::runtime_error);
  EXPECT_THROW(v.string(), std::runtime_error);
  EXPECT_THROW(v.pointer(), std::runtime_error);
}

TEST(FunctionTest, VariantFromDouble) {
  Variant v(Variant::FromDouble(42));
  EXPECT_EQ(fun::DOUBLE, v.type());
  EXPECT_EQ(42, v.number());
  EXPECT_THROW(v.string(), std::runtime_error);
  EXPECT_THROW(v.pointer(), std::runtime_error);
}

TEST(FunctionTest, VariantFromString) {
  Variant v(Variant::FromString("test"));
  EXPECT_EQ(fun::STRING, v.type());
  EXPECT_STREQ("test", v.string());
  EXPECT_THROW(v.number(), std::runtime_error);
  EXPECT_THROW(v.pointer(), std::runtime_error);
}

TEST(FunctionTest, VariantFromPointer) {
  int dummy;
  Variant v(Variant::FromPointer(&dummy));
  EXPECT_EQ(fun::POINTER, v.type());
  EXPECT_EQ(&dummy, v.pointer());
  EXPECT_THROW(v.number(), std::runtime_error);
  EXPECT_THROW(v.string(), std::runtime_error);
}

TEST(FunctionTest, MakeVariant) {
  EXPECT_EQ(42, MakeVariant(42).number());
  EXPECT_STREQ("test", MakeVariant("test").string());
  int dummy = 0;
  EXPECT_EQ(&dummy, MakeVariant(&dummy).pointer());
}

TEST(FunctionTest, CompareVariants) {
  EXPECT_TRUE(Variant() == Variant());
  EXPECT_FALSE(Variant() != Variant());
  EXPECT_TRUE(Variant::FromDouble(42) == Variant::FromDouble(42));
  EXPECT_FALSE(Variant::FromDouble(42) != Variant::FromDouble(42));

  char s1[] = "abc", s2[] = "abc", s3[] = "def";
  EXPECT_TRUE(Variant::FromString(s1) == Variant::FromString(s2));
  EXPECT_FALSE(Variant::FromString(s1) != Variant::FromString(s2));
  EXPECT_FALSE(Variant::FromString(s1) == Variant::FromString(s3));
  EXPECT_TRUE(Variant::FromString(s1) != Variant::FromString(s3));
  EXPECT_FALSE(Variant::FromString(s1) == Variant());
  EXPECT_FALSE(Variant::FromString("42") == Variant::FromDouble(42));

  EXPECT_TRUE(Variant::FromPointer(s1) == Variant::FromPointer(s1));
  EXPECT_FALSE(Variant::FromPointer(s1) != Variant::FromPointer(s1));
  EXPECT_FALSE(Variant::FromPointer(s1) == Variant::FromPointer(s2));
  EXPECT_TRUE(Variant::FromPointer(s1) != Variant::FromPointer(s2));
  EXPECT_FALSE(Variant::FromPointer(s1) == Variant());
  EXPECT_FALSE(Variant::FromString(s1) == Variant::FromPointer(s1));
}

TEST(FunctionTest, VariantOutput) {
  std::stringstream ss;
  ss << Variant::FromDouble(42);
  EXPECT_EQ("42", ss.str());
  ss.str("");
  ss << Variant::FromString("test");
  EXPECT_EQ("\"test\"", ss.str());
}

TEST(FunctionTest, EmptyTable) {
  Table t("", 0);
  EXPECT_STREQ("", t.name());
  EXPECT_EQ(0u, t.num_rows());
  EXPECT_EQ(0u, t.num_cols());
  EXPECT_THROW(t.GetColName(0), std::invalid_argument);
  EXPECT_THROW(t(0, 0), std::invalid_argument);
}

TEST(FunctionTest, TableArity) {
  Table t1("", 0);
  EXPECT_EQ(0u, t1.arity());
  Table t2("", 2);
  EXPECT_EQ(1u, t2.arity());
  Table t3("", 3);
  EXPECT_EQ(1u, t3.arity());
  t3.SetArity(3);
  EXPECT_EQ(3u, t3.arity());
  EXPECT_THROW(t3.SetArity(4), std::invalid_argument);
}

TEST(FunctionTest, Table) {
  Table t("Test", 3);
  EXPECT_STREQ("Test", t.name());
  EXPECT_EQ(0u, t.num_rows());
  EXPECT_EQ(3u, t.num_cols());
  EXPECT_THROW(t.GetColName(0), std::invalid_argument);
  EXPECT_THROW(t(0, 0), std::invalid_argument);
  t = "c1", "c2", "c3",
       11,  "v12", 13,
      "v21", 22,  "v23";
  EXPECT_THROW(t.GetColName(-1), std::invalid_argument);
  EXPECT_STREQ("c1", t.GetColName(0));
  EXPECT_STREQ("c2", t.GetColName(1));
  EXPECT_STREQ("c3", t.GetColName(2));
  EXPECT_THROW(t.GetColName(3), std::invalid_argument);
  EXPECT_EQ(11, t(0, 0).number());
  EXPECT_STREQ("v12", t(0, 1).string());
  EXPECT_EQ(13, t(0, 2).number());
  EXPECT_STREQ("v21", t(1, 0).string());
  EXPECT_EQ(22, t(1, 1).number());
  EXPECT_STREQ("v23", t(1, 2).string());
}

TEST(FunctionTest, TableComparison) {
  Table t1("t1", 0), t2("t2", 0);
  EXPECT_TRUE(t1 == t2);
  Table t3("", 2);
  EXPECT_FALSE(t1 == t3);
  Table t4("", 2);
  EXPECT_TRUE(t3 == t4);
  t3 = "a", "b",
       "c", "d";
  EXPECT_FALSE(t3 == t4);
  t4 = "a", "b",
       "c", "d";
  EXPECT_TRUE(t3 == t4);
  Table t5("", 1);
  t5 = "a", "b", "c", "d";
  EXPECT_FALSE(t3 == t5);
  Table t6("", 2);
  t6 = "a", "b";
  EXPECT_FALSE(t3 == t6);
  Table t7("", 2);
  EXPECT_FALSE(t7 == t6);
}

TEST(FunctionTest, TableOutput) {
  std::stringstream ss;
  Table t("", 2);
  t = "a", "b",
      "c", "d";
  ss << t;
  EXPECT_EQ("a b \n\"c\" \"d\" \n", ss.str());
}

TEST(FunctionTest, Library) {
  Library lib("testlib.dll");
  EXPECT_EQ(0u, lib.GetNumFunctions());
  EXPECT_TRUE(lib.GetFunction("foo") == nullptr);
  lib.Load();
  EXPECT_EQ(2u, lib.GetNumFunctions());
  EXPECT_TRUE(lib.GetFunction("nonexistent") == nullptr);
  const func_info *fi = lib.GetFunction("foo");
  EXPECT_TRUE(fi != nullptr);
  EXPECT_STREQ(fi->name, "foo");
  EXPECT_EQ(42, fi->funcp(nullptr));
  EXPECT_EQ("duplicate function 'bar'", lib.error());
  EXPECT_TRUE(lib.GetHandler("nonexistent") == nullptr);
  const Handler *handler = lib.GetHandler("testhandler");
  EXPECT_TRUE(handler != nullptr);
  Table t("", 0);
  handler->Read("", &t);
}

int LookupTest(AmplExports *, TableInfo *ti) {
  char s[] = "v2";
  char *svals[] = {0, s, 0, 0};
  double dvals[] = {42, 0, 4, 2};
  return ti->Lookup(dvals, svals, ti);
}

TEST(FunctionTest, TableLookup) {
  Library lib("");
  Handler handler(&lib, 0, LookupTest);
  Table t("", 4);
  t = "c1", "c2", "c3", "c4",
       42,  "v0",  1,    1,
       11,  "v1",  2,    1,
       42,  "v2",  3,    1,
       42,  "v2",  4,    1;
  EXPECT_EQ(0, handler.Write("", t, Handler::NOTHROW));
  t.SetArity(2);
  EXPECT_EQ(2, handler.Write("", t, Handler::NOTHROW));
  t.SetArity(3);
  EXPECT_EQ(3, handler.Write("", t, Handler::NOTHROW));
  t.SetArity(4);
  EXPECT_EQ(-1, handler.Write("", t, Handler::NOTHROW));
}

const double ITEMS[] = {5, 7, 11, 13, 17, 19, 23, 29, 31};

void CheckTuple(unsigned size, const Tuple &t) {
  EXPECT_EQ(size, t.size());
  Tuple copy(t);
  for (unsigned i = 0; i < size; ++i) {
    EXPECT_EQ(ITEMS[i], t[i].number());
    EXPECT_EQ(ITEMS[i], copy[i].number());
    copy[i] = Variant::FromDouble(42);
    EXPECT_EQ(42, copy[i].number());
  }
}

TEST(FunctionTest, MakeArgs) {
  CheckTuple(1, MakeArgs(5));
  CheckTuple(2, MakeArgs(5, 7));
  CheckTuple(3, MakeArgs(5, 7, 11));
  CheckTuple(4, MakeArgs(5, 7, 11, 13));
  CheckTuple(5, MakeArgs(5, 7, 11, 13, 17));
  CheckTuple(6, MakeArgs(5, 7, 11, 13, 17, 19));
  CheckTuple(9, MakeArgs(5, 7, 11, 13, 17, 19, 23, 29, 31));
}

TEST(FunctionTest, TupleOutput) {
  std::ostringstream os;
  os << Tuple(1, Variant::FromDouble(42));
  EXPECT_EQ("(42)", os.str());
  os.str("");
  os << MakeArgs(3, 5, 7);
  EXPECT_EQ("(3, 5, 7)", os.str());
}

void CheckBitSet(const char *expected, const BitSet &bs) {
  size_t size = std::strlen(expected);
  EXPECT_EQ(size, bs.size());
  BitSet copy(bs);
  for (size_t i = 0; i < size; ++i) {
    bool value = expected[i] != '0';
    EXPECT_EQ(value, bs[i]);
    EXPECT_EQ(value, copy[i]);
    copy[i] = !value;
    EXPECT_EQ(!value, copy[i]);
  }
}

TEST(FunctionTest, BitSet) {
  EXPECT_EQ(0u, BitSet().size());
  CheckBitSet("0", BitSet(1, false));
  CheckBitSet("11", BitSet(2, true));
  CheckBitSet("000", BitSet(3, false));
  CheckBitSet("1010", BitSet("1010"));
}

TEST(FunctionTest, FunctionWithTypes) {
  FunctionWithTypes<double> f1;
  EXPECT_EQ(1u, f1.GetNumArgs());
  EXPECT_EQ(fun::DOUBLE, f1.GetArgType(0));
  EXPECT_THROW(f1.GetArgType(1), std::out_of_range);

  FunctionWithTypes<int, double> f2;
  EXPECT_EQ(2u, f2.GetNumArgs());
  EXPECT_EQ(fun::INT, f2.GetArgType(0));
  EXPECT_EQ(fun::DOUBLE, f2.GetArgType(1));
  EXPECT_THROW(f2.GetArgType(2), std::out_of_range);

  FunctionWithTypes<double, int, double> f3;
  EXPECT_EQ(3u, f3.GetNumArgs());
  EXPECT_EQ(fun::DOUBLE, f3.GetArgType(0));
  EXPECT_EQ(fun::INT, f3.GetArgType(1));
  EXPECT_EQ(fun::DOUBLE, f3.GetArgType(2));
  EXPECT_THROW(f3.GetArgType(3), std::out_of_range);

  FunctionWithTypes<int, double, int, double> f4;
  EXPECT_EQ(4u, f4.GetNumArgs());
  EXPECT_EQ(fun::INT, f4.GetArgType(0));
  EXPECT_EQ(fun::DOUBLE, f4.GetArgType(1));
  EXPECT_EQ(fun::INT, f4.GetArgType(2));
  EXPECT_EQ(fun::DOUBLE, f4.GetArgType(3));
  EXPECT_THROW(f4.GetArgType(4), std::out_of_range);
}

double Poly1(int a) {
  return a;
}

TEST(FunctionTest, FunctionPointer1) {
  typedef fun::FunctionPointer1<int, double> Fun;
  Fun f(Poly1);
  EXPECT_EQ(1u, f.GetNumArgs());
  EXPECT_EQ(fun::INT, f.GetArgType(0));
  EXPECT_THROW(f.GetArgType(1), std::out_of_range);
  EXPECT_EQ(3, f(MakeArgs(3)));
  EXPECT_THROW(f(MakeArgs(3, 5)), std::invalid_argument);
  f = FunctionPointer(Poly1);
}

double Poly2(double a, int b) {
  return a * 10 + b;
}

TEST(FunctionTest, FunctionPointer2) {
  typedef fun::FunctionPointer2<double, int, double> Fun;
  Fun f(Poly2);
  EXPECT_EQ(2u, f.GetNumArgs());
  EXPECT_EQ(fun::DOUBLE, f.GetArgType(0));
  EXPECT_EQ(fun::INT, f.GetArgType(1));
  EXPECT_THROW(f.GetArgType(2), std::out_of_range);
  EXPECT_THROW(f(MakeArgs(3)), std::invalid_argument);
  EXPECT_EQ(35, f(MakeArgs(3, 5)));
  EXPECT_THROW(f(MakeArgs(3, 5, 7)), std::invalid_argument);
  f = FunctionPointer(Poly2);
}

double Poly3(int a, double b, unsigned c) {
  return a * 100 + b * 10 + c;
}

TEST(FunctionTest, FunctionPointer3) {
  typedef fun::FunctionPointer3<int, double, unsigned, double> Fun;
  Fun f(Poly3);
  EXPECT_EQ(3u, f.GetNumArgs());
  EXPECT_EQ(fun::INT, f.GetArgType(0));
  EXPECT_EQ(fun::DOUBLE, f.GetArgType(1));
  EXPECT_EQ(fun::UINT, f.GetArgType(2));
  EXPECT_THROW(f.GetArgType(3), std::out_of_range);
  EXPECT_THROW(f(MakeArgs(3, 5)), std::invalid_argument);
  EXPECT_EQ(357, f(MakeArgs(3, 5, 7)));
  EXPECT_THROW(f(MakeArgs(3, 5, 7, 11)), std::invalid_argument);
  f = FunctionPointer(Poly3);
}

double Poly4(double a, int b, double c, unsigned d) {
  return a * 1000 + b * 100 + c * 10 + d;
}

TEST(FunctionTest, FunctionPointer4) {
  typedef fun::FunctionPointer4<double, int, double, unsigned, double> Fun;
  Fun f(Poly4);
  EXPECT_EQ(4u, f.GetNumArgs());
  EXPECT_EQ(fun::DOUBLE, f.GetArgType(0));
  EXPECT_EQ(fun::INT, f.GetArgType(1));
  EXPECT_EQ(fun::DOUBLE, f.GetArgType(2));
  EXPECT_EQ(fun::UINT, f.GetArgType(3));
  EXPECT_THROW(f.GetArgType(4), std::out_of_range);
  EXPECT_THROW(f(MakeArgs(2, 3, 5)), std::invalid_argument);
  EXPECT_EQ(2357, f(MakeArgs(2, 3, 5, 7)));
  EXPECT_THROW(f(MakeArgs(2, 3, 5, 7, 11)), std::invalid_argument);
  f = FunctionPointer(Poly4);
}

double Poly5(double a, int b, unsigned c, double d, int e) {
  return a * 10000 + b * 1000 + c * 100 + d * 10 + e;
}
TEST(FunctionTest, FunctionPointer5) {
  typedef fun::FunctionPointer5<
      double, int, unsigned, double, int, double> Fun;
  Fun f(Poly5);
  EXPECT_EQ(5u, f.GetNumArgs());
  EXPECT_EQ(fun::DOUBLE, f.GetArgType(0));
  EXPECT_EQ(fun::INT, f.GetArgType(1));
  EXPECT_EQ(fun::UINT, f.GetArgType(2));
  EXPECT_EQ(fun::DOUBLE, f.GetArgType(3));
  EXPECT_EQ(fun::INT, f.GetArgType(4));
  EXPECT_THROW(f.GetArgType(5), std::out_of_range);
  EXPECT_THROW(f(MakeArgs(1, 2, 3, 5)), std::invalid_argument);
  EXPECT_EQ(12357, f(MakeArgs(1, 2, 3, 5, 7)));
  EXPECT_THROW(f(MakeArgs(1, 2, 3, 5, 7, 11)), std::invalid_argument);
  f = FunctionPointer(Poly5);
}

TEST(FunctionTest, BindOne) {
  typedef fun::FunctionPointer3<int, double, unsigned, double> Fun;
  fun::OneBinder<Fun, double> f(Fun(Poly3), 5, 1);
  EXPECT_EQ(2u, f.GetNumArgs());
  EXPECT_EQ(fun::INT, f.GetArgType(0));
  EXPECT_EQ(fun::UINT, f.GetArgType(1));
  EXPECT_THROW(f.GetArgType(2), std::out_of_range);
  EXPECT_EQ(357, f(MakeArgs(3, 7)));
  f = fun::BindOne(Fun(Poly3), 5.0, 1);
}

#ifdef HAVE_LAMBDAS
TEST(FunctionTest, Fun) {
  double (*pow)(double, double) = &std::pow;
  Differentiator dx, dy;
  EXPECT_NEAR(4.77259,
      dx([&](double x) { return dy(bind1st(ptr_fun(pow), x), 2); }, 2),
      1e-5);
}
#endif

double TestFun2(const Tuple &args) {
  return args[0].number() * 10 + args[1].number();
}

double TestFun3(const Tuple &args) {
  return args[0].number() * 100 + args[1].number() * 10 + args[2].number();
}

TEST(FunctionTest, BindAllButOne) {
  typedef double (*Fun)(const Tuple &);

  EXPECT_EQ(35, BindAllButOne(TestFun2, MakeArgs(0, 5), 0)(3));
  EXPECT_EQ(35, fun::AllButOneBinder<Fun>(TestFun2, MakeArgs(0, 5), 0)(3));
  EXPECT_EQ(53, BindAllButOne(TestFun2, MakeArgs(5, 0), 1)(3));
  EXPECT_EQ(53, fun::AllButOneBinder<Fun>(TestFun2, MakeArgs(5, 0), 1)(3));

  EXPECT_THROW(
      fun::AllButOneBinder<Fun>(TestFun2, MakeArgs(0, 5), 2),
      std::out_of_range);

  EXPECT_EQ(357, BindAllButOne(TestFun3, MakeArgs(0, 5, 7), 0)(3));
  EXPECT_EQ(357, fun::AllButOneBinder<Fun>(TestFun3, MakeArgs(0, 5, 7), 0)(3));
  EXPECT_EQ(537, BindAllButOne(TestFun3, MakeArgs(5, 0, 7), 1)(3));
  EXPECT_EQ(537, fun::AllButOneBinder<Fun>(TestFun3, MakeArgs(5, 0, 7), 1)(3));
  EXPECT_EQ(573, BindAllButOne(TestFun3, MakeArgs(5, 7, 0), 2)(3));
  EXPECT_EQ(573, fun::AllButOneBinder<Fun>(TestFun3, MakeArgs(5, 7, 0), 2)(3));

  EXPECT_THROW(
      fun::AllButOneBinder<Fun>(TestFun3, MakeArgs(0, 5, 7), 3),
      std::out_of_range);
}

double Hypot(double x, double y) {
  return sqrt(x * x + y * y);
}

typedef double (*DoubleFun)(double x);
DoubleFun GetDoubleFun(DoubleFun f) { return f; }

TEST(FunctionTest, Differentiator) {
  Differentiator diff;
  double error = std::numeric_limits<double>::quiet_NaN();
  EXPECT_NEAR(1, diff(GetDoubleFun(std::sin), 0, &error), 1e-7);
  EXPECT_NEAR(0, error, 1e-10);
  EXPECT_NEAR(0.25, diff(GetDoubleFun(sqrt), 4), 1e-7);
  EXPECT_NEAR(0, diff(std::bind2nd(ptr_fun(Hypot), -5), 0), 1e-7);
}

TEST(FunctionTest, DifferentiatorPropagatesNaN) {
  Differentiator diff;
  EXPECT_NE(0, isnan(sqrt(-1.0)));
  EXPECT_NE(0, isnan(diff(GetDoubleFun(sqrt), -1)));
}

TEST(FunctionTest, DifferentiatorDetectsNaN) {
  Differentiator diff;
  EXPECT_EQ(0, std::bind2nd(ptr_fun(Hypot), 0)(0));
  EXPECT_NE(0, isnan(diff(std::bind2nd(ptr_fun(Hypot), 0), 0)));
  EXPECT_EQ(-std::numeric_limits<double>::infinity(), std::log(0.0));
  EXPECT_NE(0, isnan(diff(GetDoubleFun(std::log), 0)));
}

double PositiveOrNaN(double x) {
  return x >= 0 ? x : std::numeric_limits<double>::quiet_NaN();
}

TEST(FunctionTest, DifferentiatorRightDeriv) {
  // Differentiator should use the right derivative if the function is not
  // defined for a negative argument.
  Differentiator diff;
  EXPECT_NE(0, isnan(PositiveOrNaN(-1e-7)));
  EXPECT_NEAR(1, diff(ptr_fun(PositiveOrNaN), 0), 1e-7);
}

TEST(FunctionTest, FunctionInfoArgNames) {
  FunctionInfo fi;
  EXPECT_THROW(fi.GetArgName(0), std::out_of_range);
  fi.SetArgNames("x y z");
  EXPECT_EQ("x", fi.GetArgName(0));
  EXPECT_EQ("y", fi.GetArgName(1));
  EXPECT_EQ("z", fi.GetArgName(2));
  EXPECT_THROW(fi.GetArgName(3), std::out_of_range);
}

struct TestFunctionInfo : FunctionInfo {
  Result GetDerivative(const Function &, unsigned, const Tuple &) const {
    return Result(42);
  }

  Result GetSecondDerivative(
      const Function &, unsigned, unsigned, const Tuple &) const {
    return Result(11);
  }
};

TEST(FunctionTest, FunctionInfoGetDerivative) {
  Library lib("");
  FunctionInfo fi1;
  Function f(&lib, 0, 0);
  EXPECT_NE(0, isnan(fi1.GetDerivative(f, 0, MakeArgs(0)).value()));
  EXPECT_NE(0, isnan(fi1.GetSecondDerivative(f, 0, 0, MakeArgs(0)).value()));
  TestFunctionInfo fi2;
  EXPECT_EQ(42, fi2.GetDerivative(f, 0, MakeArgs(0)).value());
  EXPECT_EQ(11, fi2.GetSecondDerivative(f, 0, 0, MakeArgs(0)).value());
}

TEST(FunctionTest, FunctionInfoResult) {
  EXPECT_EQ(42, FunctionInfo::Result(42).value());
  EXPECT_TRUE(FunctionInfo::Result().error() == nullptr);
  EXPECT_NE(0, isnan(FunctionInfo::Result().value()));
  EXPECT_NE(0, isnan(FunctionInfo::Result("oops").value()));
  EXPECT_STREQ("oops", FunctionInfo::Result("oops").error());
}

TEST(FunctionTest, FunctionResult) {
  static const real ARGS[] = {5, 7, 11, 13, 17};
  Function::Result r(42, vector<real>(ARGS, ARGS + 2),
      vector<real>(ARGS + 2, ARGS + 5), nullptr);
  EXPECT_EQ(42, r);
  EXPECT_EQ(5, r.deriv());
  EXPECT_EQ(5, r.deriv(0));
  EXPECT_EQ(7, r.deriv(1));
  EXPECT_THROW(r.deriv(2), std::out_of_range);
  EXPECT_EQ(11, r.hes());
  EXPECT_EQ(11, r.hes(0));
  EXPECT_EQ(13, r.hes(1));
  EXPECT_EQ(17, r.hes(2));
  EXPECT_THROW(r.hes(3), std::out_of_range);
  EXPECT_TRUE(r.error() == nullptr);
}

TEST(FunctionTest, FunctionResultError) {
  static const real ARGS[] = {5, 7, 11, 13, 17};
  const char *error = "brain overflow";
  Function::Result r(42, vector<real>(ARGS, ARGS + 2),
      vector<real>(ARGS + 2, ARGS + 5), error);
  EXPECT_THROW(static_cast<double>(r), std::runtime_error);
  EXPECT_THROW(r.deriv(), std::runtime_error);
  EXPECT_THROW(r.deriv(0), std::runtime_error);
  EXPECT_THROW(r.hes(), std::runtime_error);
  EXPECT_THROW(r.hes(0), std::runtime_error);
  EXPECT_STREQ(error, r.error());
}

struct CallData {
  AmplExports *ae;
  TMInfo *tmi;
  int n;
  int nr;
  vector<real> ra;
  real *derivs;
  real *hes;
  char *dig;
  char *error;
};

real Test(arglist *args) {
  CallData *data = reinterpret_cast<CallData*>(args->funcinfo);
  data->ae = args->AE;
  data->tmi = args->TMI;
  data->n = args->n;
  data->nr = args->nr;
  data->ra = vector<real>(args->ra, args->ra + args->n);
  data->derivs = args->derivs;
  data->hes = args->hes;
  data->dig = args->dig;
  data->error = args->Errmsg;
  if (args->ra[0] < 0)
    args->Errmsg = const_cast<char*>("oops");
  if (args->derivs) {
    args->derivs[0] = 123;
    args->derivs[1] = 456;
    if (args->n > 2)
      args->derivs[2] = 789;
  }
  if (args->hes) {
    args->hes[0] = 12;
    args->hes[1] = 34;
    args->hes[2] = 56;
  }
  return 42;
}

class TestFunction {
 private:
  Library lib_;
  func_info fi_;
  Function f_;

 public:
  explicit TestFunction(int nargs)
  : lib_(""), fi_(), f_(&lib_, &fi_, 0) {
    fi_.nargs = nargs;
    fi_.funcp = Test;
  }

  const Function &get() const { return f_; }
};

TEST(FunctionTest, FunctionCall) {
  TestFunction f(1);
  CallData data = {};
  EXPECT_EQ(42, f.get()(777, 0, BitSet(), &data));
  ASSERT_EQ(1, data.n);
  EXPECT_EQ(1, data.nr);
  EXPECT_EQ(777, data.ra[0]);
  EXPECT_TRUE(data.derivs == nullptr);
  EXPECT_TRUE(data.hes == nullptr);
  EXPECT_TRUE(data.dig == nullptr);
  EXPECT_TRUE(data.error == nullptr);
}

TEST(FunctionTest, FunctionReturnsError) {
  TestFunction f(1);
  CallData data = {};
  EXPECT_STREQ("oops", f.get()(-1, 0, BitSet(), &data).error());
}

TEST(FunctionTest, FunctionReturnsDerivs) {
  TestFunction f(3);
  CallData data = {};
  Function::Result res =
      f.get()(MakeArgs(11, 22, 33), fun::DERIVS, BitSet(), &data);
  EXPECT_EQ(42, res);
  ASSERT_EQ(3, data.n);
  EXPECT_EQ(3, data.nr);
  EXPECT_EQ(11, data.ra[0]);
  EXPECT_EQ(22, data.ra[1]);
  EXPECT_EQ(33, data.ra[2]);
  EXPECT_EQ(123, res.deriv(0));
  EXPECT_EQ(456, res.deriv(1));
  EXPECT_EQ(789, res.deriv(2));
  EXPECT_THROW(res.deriv(3), std::out_of_range);
  EXPECT_TRUE(data.hes == nullptr);
  EXPECT_TRUE(data.dig == nullptr);
  EXPECT_TRUE(data.error == nullptr);
}

TEST(FunctionTest, FunctionReturnsHes) {
  TestFunction f(2);
  CallData data = {};
  Function::Result res = f.get()(MakeArgs(111, 222), fun::HES, BitSet(), &data);
  EXPECT_EQ(42, res);
  ASSERT_EQ(2, data.n);
  EXPECT_EQ(2, data.nr);
  EXPECT_EQ(111, data.ra[0]);
  EXPECT_EQ(222, data.ra[1]);
  EXPECT_EQ(123, res.deriv(0));
  EXPECT_EQ(456, res.deriv(1));
  EXPECT_THROW(res.deriv(2), std::out_of_range);
  EXPECT_EQ(12, res.hes(0));
  EXPECT_EQ(34, res.hes(1));
  EXPECT_EQ(56, res.hes(2));
  EXPECT_THROW(res.hes(3), std::out_of_range);
  EXPECT_TRUE(data.dig == nullptr);
  EXPECT_TRUE(data.error == nullptr);
}

TEST(FunctionTest, FunctionArgNames) {
  Library lib("");
  FunctionInfo fi;
  Function f(&lib, 0, &fi);
  EXPECT_THROW(f.GetArgName(0), std::out_of_range);
  fi.SetArgNames("x y z");
  EXPECT_EQ("x", f.GetArgName(0));
  EXPECT_EQ("y", f.GetArgName(1));
  EXPECT_EQ("z", f.GetArgName(2));
  EXPECT_THROW(f.GetArgName(3), std::out_of_range);
}

TEST(FunctionTest, FunctifonGetDerivative) {
  Library lib("");
  FunctionInfo fi1;
  Function f1(&lib, 0, &fi1);
  EXPECT_NE(0, isnan(f1.GetDerivative(0, MakeArgs(0)).value()));
  EXPECT_NE(0, isnan(f1.GetSecondDerivative(0, 0, MakeArgs(0)).value()));
  TestFunctionInfo fi2;
  Function f2(&lib, 0, &fi2);
  EXPECT_EQ(42, f2.GetDerivative(0, MakeArgs(0)).value());
  EXPECT_EQ(11, f2.GetSecondDerivative(0, 0, MakeArgs(0)).value());
}

real ASLHypot(arglist *al) {
  double x = al->ra[0];
  double y = al->ra[1];
  real result = sqrt(x * x + y * y);
  if (al->derivs) {
    real *derivs = al->derivs;
    derivs[0] = x / result;
    derivs[1] = y / result;
  }
  return result;
}

TEST(FunctionTest, DerivativeBinder) {
  Library lib("");
  func_info fi = {};
  fi.nargs = 2;
  fi.funcp = ASLHypot;
  Function f(&lib, &fi, 0);
  DerivativeBinder d(f, 0, 1, MakeArgs(1, 0));
  ASSERT_EQ(1, d(0));
  ASSERT_NEAR(1 / sqrt(2.0), d(1), 1e-10);
  d = DerivativeBinder(f, 1, 1, MakeArgs(1, 0));
  ASSERT_EQ(0, d(0));
  ASSERT_NEAR(1 / sqrt(2.0), d(1), 1e-10);
  EXPECT_THROW(DerivativeBinder(f, 2, 0, MakeArgs(0, 0)), std::out_of_range);
  EXPECT_THROW(DerivativeBinder(f, 0, 2, MakeArgs(0, 0)), std::out_of_range);
}
}
