/**
 * Example using "easy" NLWriter2 API, mp::NLSOL_Easy,
 * to solve MIQP.
 *
 * For full NLWriter2 API, see mp::NLSOL
 * and related tests/examples.
 */

#include <cstdlib>
#include <cstring>
#include <cmath>

#include "mp/nlsol-easy.h"

/// Helper class to build model
/// and check solution
class ModelBuilder {
public:
  /// Obtain model.
  /// Valid with *this's lifetime.
  mp::NLModel_Easy GetModel() {
    mp::NLModel_Easy nlme(prob_name_);

    nlme.SetCols({(int)var_lb_.size(),
                  var_lb_.data(), var_ub_.data(),
                  var_type_.data()});
    nlme.SetColNames(var_names_.data());

    nlme.SetRows(row_lb_.size(), row_lb_.data(), row_ub_.data(),
                 {(int)row_lb_.size(), A_format_,
                  A_index_.size(),
                  A_start_.data(), A_index_.data(), A_value_.data()});
    nlme.SetRowNames(row_names_.data());

    nlme.SetLinearObjective(obj_sense_, obj_c0_, obj_c_.data());
    nlme.SetHessian(Q_format_,
                    {(int)var_lb_.size(), 0,
                     Q_index_.size(),
                     Q_start_.data(), Q_index_.data(), Q_value_.data()});
    nlme.SetObjName(obj_name_);

    return nlme;
  }

  /// Check solution
  bool Check(mp::NLSOL_Easy::Solution sol) {
    if (!ApproxEqual(sol.obj_val_, obj_val_ref_)) {
      printf("MIQP 1: wrong obj val (%.17g !~ %.17g)\n",
             sol.obj_val_, obj_val_ref_);
      return false;
    }
    for (auto i=sol.x_.size(); i--; )
      if (!ApproxEqual(x_ref_[i], sol.x_[i])) {
        printf("MIQP 1: wrong x[%ld] (%.17g !~ %.17g)\n",
               i+1, sol.x_[i], x_ref_[i]);
        return false;
      }
    printf("MIQP 1: solution check ok.\n");
    return true;
  }

protected:
  static bool ApproxEqual(double n, double m) {
    return std::fabs(n-m)
        <= 1e-5 * std::min(1.0, std::fabs(n)+std::fabs(m));
  }

private:
  /*
var x1_4 >=0 <=0;
var x2_5 >=-3 <=20 integer;
var x3_6 binary;
var x4_3 >=-1 <=1e20 integer;
var x5_1 >=-1 <=-1;
var x6_2 >=-2 <=10;

minimize Z:
   x2_5 + 3.24
     + 5*x4_3*x4_3 + 6*x4_3*x6_2 + 7*x5_1^2; ## + x6^2?

s.t. C1: x2_5 + x3_6 + x4_3 + x6_2 == 15;

s.t. C2: x2_5 -x3_6 -x4_3 + x6_2 >= 10;
   */
  const char* prob_name_ {"NLSOL_Easy_model"};
  std::vector<double> var_lb_
  {0, -3, 0, -1, -1, -2};
  std::vector<double> var_ub_
  {0, 20, 1, 1e20, -1, 10};
  std::vector<int> var_type_
  {0, 1, 1, 1, 0, 0};
  std::vector<const char*> var_names_
  {"x1_4", "x2_6", "x3_5", "x4_3", "x5_1", "x6_2"};
  int A_format_ {NLW2_MatrixFormatRowwise};
  std::vector<size_t> A_start_ {0, 4};
  std::vector<int> A_index_ {1,2,3,5,1,2,3,5};
  std::vector<double> A_value_{1,1,1,1,1,-1,-1,1};
  std::vector<double> row_lb_{15, 10};
  std::vector<double> row_ub_{15, INFINITY};
  std::vector<const char*> row_names_ {"C1", "C2"};
  int obj_sense_ {0};
  double obj_c0_ {3.24};
  std::vector<double> obj_c_ {0,1,0,0,0,0};
  int Q_format_ {NLW2_HessianFormatSquare};
  std::vector<size_t> Q_start_ {0,0,0,0,2,3};
  std::vector<int> Q_index_ {3, 5, 4};
  std::vector<double> Q_value_ {10, 12, 14};
  const char* obj_name_ {"obj[1]"};

  /// Solution
  std::vector<double> x_ref_
  {0, 5, 1, -1, -1, 10};
  double obj_val_ref_ {-39.76};
};

/// Solver with given parameters
bool SolveAndCheck(std::string solver, std::string sopts,
                   bool binary, std::string stub) {
  ModelBuilder mdlbld;
  auto nlme = mdlbld.GetModel();
  mp::NLSOL_Easy nlse;
  auto nlopts = NLW2_Make_NLOptionsBasic_C_Default();
  nlopts.n_text_mode_ = !binary;
  nlopts.want_nl_comments_ = 1;
  nlse.SetNLOptions(nlopts);
  nlse.SetFileStub(stub);
  if (auto sol = nlse.Solve(nlme, solver, sopts)) {
    if (!mdlbld.Check(sol)) {
      printf("Solution check failed.\n");
      return false;
    }
  } else {
    printf("%s\n", nlse.GetErrorMessage());
    return false;
  }

  return true;
}

/// Invoke:
///   (this_exe) ipopt ["outlev=1 timelim=500" [text [filestub]]]
int main(int argc, const char* const* argv) {
  if (argc<2) {
    printf("%s\n",
           "AMPL NL writer example.\n"
           "Usage:\n"
           "  <this_exe> <solver> [\"<solver_options>\" [binary/text [<stub>]]],\n\n"
           "where <solver> is ipopt, gurobi, minos, ...;\n"
           "binary/text is the NL format (default: binary.)\n"
           "Examples:\n"
           "  <this_exe> ipopt \"\" text /tmp/stub\n"
           "  <this_exe> gurobi \"nonconvex=2 funcpieces=-2 funcpieceratio=1e-4\"");
    exit(0);
  }
  std::string solver = (argc>1) ? argv[1] : "minos";
  std::string sopts = (argc>2) ? argv[2] : "";
  bool binary = (argc<=3) || std::strcmp("text", argv[3]);
  std::string stub = (argc>4) ? argv[4] : "";

  if (!SolveAndCheck(solver, sopts, binary, stub))
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}
