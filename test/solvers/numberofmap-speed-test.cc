#include "ilogcp/concert.h"
#include "mp/utils-clock.h"
#include "mp/problem.h"

struct CreateVar {
  int operator()() { return 0; }
};

int main() {
  mp::Problem p;
  int num_exprs = 10000;
  mp::NumberOfMap<int, CreateVar> map((CreateVar()));
  mp::NumericConstant n = p.MakeNumericConstant(0);
  std::vector<mp::IteratedExpr> exprs(num_exprs);
  for (int i = 0; i < num_exprs; ++i) {
    mp::Problem::NumberOfExprBuilder b = p.BeginNumberOf(2, n);
    b.AddArg(p.MakeVariable(i));
    exprs[i] = p.EndNumberOf(b);
  }
  std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  for (int i = 0; i < num_exprs; ++i)
    map.Add(0, exprs[i]);
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  fmt::print("Executed NumberOfMap.Add {} times in {} s.\n",
    num_exprs, mp::duration_cast< mp::duration<double> >(end - start).count());
}
