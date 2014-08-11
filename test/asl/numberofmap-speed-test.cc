#include "solvers/util/aslbuilder.h"
#include "solvers/util/clock.h"
#include "solvers/util/nl.h"

struct CreateVar {
  int operator()() { return 0; }
};

int main() {
  ampl::internal::ASLBuilder b;
  int num_exprs = 10000;
  ampl::NLHeader h = {};
  h.num_vars = num_exprs;
  h.num_objs = 1;
  b.BeginBuild("", h, 0);
  ampl::NumberOfMap<int, CreateVar> map((CreateVar()));
  ampl::NumericConstant n = b.MakeNumericConstant(0);
  std::vector<ampl::NumberOfExpr> exprs(num_exprs);
  for (int i = 0; i < num_exprs; ++i) {
    ampl::NumericExpr args[] = {n, b.MakeVariable(i)};
    exprs[i] = b.MakeNumberOf(args);
  }
  ampl::steady_clock::time_point start = ampl::steady_clock::now();
  for (int i = 0; i < num_exprs; ++i)
    map.Add(0, exprs[i]);
  ampl::steady_clock::time_point end = ampl::steady_clock::now();
  fmt::print("Executed NumberOfMap.Add {} times in {} s.\n",
    num_exprs, ampl::duration_cast< ampl::duration<double> >(end - start).count());
}
