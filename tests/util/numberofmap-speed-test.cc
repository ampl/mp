#include "solvers/util/clock.h"
#include "tests/expr-builder.h"

struct CreateVar {
  int operator()() { return 0; }
};

int main() {
  ampl::ExprBuilder eb;
  ampl::NumberOfMap<int, CreateVar> map((CreateVar()));
  ampl::NumericConstant n = eb.AddNum(0);
  int num_exprs = 10000;
  std::vector<ampl::NumberOfExpr> exprs(num_exprs);
  for (int i = 0; i < num_exprs; ++i)
    exprs[i] = eb.AddNumberOf(n, eb.AddVar(i));
  ampl::steady_clock::time_point start = ampl::steady_clock::now();
  for (int i = 0; i < num_exprs; ++i)
    map.Add(0, exprs[i]);
  ampl::steady_clock::time_point end = ampl::steady_clock::now();
  fmt::print("Executed NumberOfMap.Add {} times in {} s.\n",
    num_exprs, ampl::duration_cast< ampl::duration<double> >(end - start).count());
}
