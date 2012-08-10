// Tests for the functional stuff.

#include <limits>
#include <sstream>
#include <cmath>
#include <cstring>

#include "gtest/gtest.h"
#include "tests/function.h"
#include "solvers/asl.h"

using std::ptr_fun;
using std::sqrt;
using std::vector;

using fun::BitSet;
using fun::Differentiator;
using fun::Function;
using fun::Tuple;

namespace {

const double ITEMS[] = {5, 7, 11, 13, 17, 19, 23, 29, 31};

void CheckTuple(unsigned size, const Tuple &t) {
  EXPECT_EQ(size, t.size());
  Tuple copy(t);
  for (unsigned i = 0; i < size; ++i) {
    EXPECT_EQ(ITEMS[i], t[i]);
    EXPECT_EQ(ITEMS[i], copy[i]);
    copy[i] = 42;
    EXPECT_EQ(42, copy[i]);
  }
}

TEST(FunctionTest, Tuple) {
  CheckTuple(1, Tuple(5));
  CheckTuple(2, Tuple(5, 7));
  CheckTuple(3, Tuple(5, 7, 11));
  CheckTuple(4, Tuple(5, 7, 11, 13));
  CheckTuple(5, Tuple(5, 7, 11, 13, 17));
  CheckTuple(6, Tuple(5, 7, 11, 13, 17, 19));
  CheckTuple(9, Tuple(5, 7, 11, 13, 17, 19, 23, 29, 31));
}

TEST(FunctionTest, TupleOutput) {
  std::ostringstream os;
  os << Tuple(42);
  EXPECT_EQ("(42)", os.str());
  os.str("");
  os << Tuple(3, 5, 7);
  EXPECT_EQ("(3, 5, 7)", os.str());
}

void CheckBitSet(const char *expected, const BitSet &bs) {
  size_t size = std::strlen(expected);
  EXPECT_EQ(size, bs.size());
  BitSet copy(bs);
  for (size_t i = 0; i < size; ++i) {
    bool value = expected[i] - '0';
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

#define EXPECT_TYPE(expected, actual) \
  { expected *e = static_cast<actual*>(0); e = e; }

TEST(FunctionTest, TernaryFunction) {
  typedef fun::ternary_function<bool, char, int, double> Fun;
  EXPECT_TYPE(bool, Fun::first_argument_type);
  EXPECT_TYPE(char, Fun::second_argument_type);
  EXPECT_TYPE(int, Fun::third_argument_type);
  EXPECT_TYPE(double, Fun::result_type);
}

double TestFun1(bool b, char c, int i) {
  EXPECT_EQ(true, b);
  EXPECT_EQ('a', c);
  EXPECT_EQ(42, i);
  return 777;
}

TEST(FunctionTest, PointerToTernaryFunction) {
  typedef fun::pointer_to_ternary_function<bool, char, int, double> Fun;
  EXPECT_TYPE(bool, Fun::first_argument_type);
  EXPECT_TYPE(char, Fun::second_argument_type);
  EXPECT_TYPE(int, Fun::third_argument_type);
  EXPECT_TYPE(double, Fun::result_type);
  Fun f(TestFun1);
  EXPECT_EQ(777, f(true, 'a', 42));
}

TEST(FunctionTest, Fun) {
  Differentiator dx, dy;
  EXPECT_NEAR(4.77259,
      dx([&](double x) { return dy(bind1st(ptr_fun(pow), x), 2); }, 2),
      1e-5);
}

double TestFun2(double x, double y, double z) {
  return x * 100 + y * 10 + z;
}

TEST(FunctionTest, Bind2Of3) {
  typedef double (*Fun)(double, double, double);
  fun::Binder2Of3<Fun> f1(TestFun2, Tuple(0, 5, 7), 0);
  EXPECT_EQ(357, Bind2Of3(TestFun2, Tuple(0, 5, 7), 0)(3));
  EXPECT_EQ(357, f1(3));
  fun::Binder2Of3<Fun> f2(TestFun2, Tuple(5, 0, 7), 1);
  EXPECT_EQ(537, Bind2Of3(TestFun2, Tuple(5, 0, 7), 1)(3));
  EXPECT_EQ(537, f2(3));
  fun::Binder2Of3<Fun> f3(TestFun2, Tuple(5, 7, 0), 2);
  EXPECT_EQ(573, Bind2Of3(TestFun2, Tuple(5, 7, 0), 2)(3));
  EXPECT_EQ(573, f3(3));

  EXPECT_THROW(fun::Binder2Of3<Fun>(TestFun2, Tuple(1), 0),
      std::out_of_range);
  EXPECT_THROW(fun::Binder2Of3<Fun>(TestFun2, Tuple(1, 2), 0),
      std::out_of_range);
  fun::Binder2Of3<Fun>(TestFun2, Tuple(1, 2, 3), 0);
  EXPECT_THROW(fun::Binder2Of3<Fun>(TestFun2, Tuple(1, 2, 3, 4), 0),
      std::out_of_range);
  EXPECT_THROW(fun::Binder2Of3<Fun>(TestFun2, Tuple(1, 2, 3), 3),
      std::out_of_range);
}

double Hypot(double x, double y) {
  return sqrt(x * x + y * y);
}

typedef double (*DoubleFun)(double);
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
  EXPECT_TRUE(std::isnan(sqrt(-1)));
  EXPECT_TRUE(std::isnan(diff(GetDoubleFun(sqrt), -1)));
}

TEST(FunctionTest, DifferentiatorDetectsNaN) {
  Differentiator diff;
  EXPECT_EQ(0, std::bind2nd(ptr_fun(Hypot), 0)(0));
  EXPECT_TRUE(std::isnan(diff(std::bind2nd(ptr_fun(Hypot), 0), 0)));
  EXPECT_EQ(-std::numeric_limits<double>::infinity(), std::log(0));
  EXPECT_TRUE(std::isnan(diff(GetDoubleFun(std::log), 0)));
}

double PositiveOrNaN(double x) {
  return x >= 0 ? x : std::numeric_limits<double>::quiet_NaN();
}

TEST(FunctionTest, DifferentiatorRightDeriv) {
  // Differentiator should use the right derivative if the function is not
  // defined for a negative argument.
  Differentiator diff;
  EXPECT_TRUE(std::isnan(PositiveOrNaN(-1e-7)));
  EXPECT_NEAR(1, diff(ptr_fun(PositiveOrNaN), 0), 1e-7);
}

TEST(FunctionTest, Result) {
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

TEST(FunctionTest, ErrorResult) {
  static const real ARGS[] = {5, 7, 11, 13, 17};
  const char *error = "brain overflow";
  Function::Result r(42, vector<real>(ARGS, ARGS + 2),
      vector<real>(ARGS + 2, ARGS + 5), error);
  EXPECT_THROW((double)r, std::runtime_error);
  EXPECT_THROW(r.deriv(), std::runtime_error);
  EXPECT_THROW(r.deriv(0), std::runtime_error);
  EXPECT_THROW(r.hes(), std::runtime_error);
  EXPECT_THROW(r.hes(0), std::runtime_error);
  EXPECT_STREQ(error, r.error());
}

struct CallData {
  AmplExports *ae;
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
  ASL testASL_;
  AmplExports ae_;
  func_info fi_;
  Function f_;

 public:
  TestFunction(int nargs)
  : testASL_(), ae_(), fi_(), f_(&testASL_, &fi_, 0) {
    testASL_.i.ae = &ae_;
    fi_.nargs = nargs;
    fi_.funcp = Test;
  }

  const AmplExports* ae() const { return &ae_; }
  const Function& get() const { return f_; }
};

TEST(FunctionTest, FunctionCall) {
  TestFunction f(1);
  CallData data = {};
  EXPECT_EQ(42, f.get()(777, 0, BitSet(), &data));
  EXPECT_EQ(f.ae(), data.ae);
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
      f.get()(Tuple(11, 22, 33), fun::DERIVS, BitSet(), &data);
  EXPECT_EQ(42, res);
  EXPECT_EQ(f.ae(), data.ae);
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
  Function::Result res = f.get()(Tuple(111, 222), fun::HES, BitSet(), &data);
  EXPECT_EQ(42, res);
  EXPECT_EQ(f.ae(), data.ae);
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

}
