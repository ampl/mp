/*
 A mathematical optimization solver.

 Copyright (C) 2012 AMPL Optimization Inc

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

extern "C" {
#include "solvers/getstub.h"
}

#undef Char
#undef filename
#undef ampl_vbtol
#undef write_sol

#include "mp/clock.h"
#include "mp/sol.h"
#include "aslsolver.h"
#include "problem.h"

using mp::SuffixView;

namespace {

// Flags for PrintUsage.
enum {
  // Print usage to stderr.
  USE_STDERR   = 1,
  // Make "-?" show some command-line options that are relevant when
  // importing functions.
  FUNC_OPTIONS = 2
};

void PrintUsage(const mp::Solver &solver, unsigned flags) {
  FILE *f = (flags & USE_STDERR) != 0 ? stderr : stdout;
  fmt::print(f, "usage: {} [options] stub [-AMPL] [<assignment> ...]\n",
             solver.name());
  fmt::print(f, "\nOptions:\n");
  static const char *options[] = {
    "-",  "end of options",
    "=",  "show name= possibilities",
    "?",  "show usage",
    "e",  "suppress echoing of assignments",
    "ix", "import user-defined functions from x; -i? gives details",
    "s",  "write .sol file (without -AMPL)",
    "u",  "just show available user-defined functions",
    "v",  "just show version",
    0
  };
  for (const char **s = options; *s; s += 2) {
    // Filter out -ix & -u options if FUNC_OPTIONS flags is not set.
    char c = **s;
    if ((flags & FUNC_OPTIONS) == 0 && (c == 'i' || c == 'u'))
      continue;
    fmt::print(f, "\t-{:3}{{{}}}\n", s[0], s[1]);
  }
}

// Adapts an ASL solution for WriteSol.
class ASLSolutionAdapter {
 private:
  const mp::Problem &problem_;
  const char *message_;
  mp::ArrayRef<fint> options_;
  mp::ArrayRef<double> values_;
  mp::ArrayRef<double> dual_values_;

 public:
  ASLSolutionAdapter(const mp::Problem &p, const char *message,
                     mp::ArrayRef<fint> options, mp::ArrayRef<double> values,
                     mp::ArrayRef<double> dual_values)
    : problem_(p), message_(message), options_(options),
      values_(values), dual_values_(dual_values) {}

  const char *message() const { return message_; }

  int num_options() const { return options_.size(); }
  int option(int index) const { return options_[index]; }

  int num_values() const { return values_.size(); }
  int value(int index) const { return values_[index]; }

  int num_dual_values() const { return dual_values_.size(); }
  int dual_value(int index) const { return dual_values_[index]; }

  SuffixView var_suffixes() const { return problem_.var_suffixes(); }
  SuffixView con_suffixes() const { return problem_.con_suffixes(); }
  SuffixView obj_suffixes() const { return problem_.obj_suffixes(); }
  SuffixView problem_suffixes() const { return problem_.problem_suffixes(); }
};
}  // namespace

void mp::SolutionWriter::HandleFeasibleSolution(
    fmt::StringRef message, const double *values,
    const double *dual_values, double) {
  ++num_solutions_;
  const char *solution_stub = solver_.solution_stub();
  if (!*solution_stub)
    return;
  fmt::Writer w;
  w << solution_stub << num_solutions_ << ".sol";
  Option_Info option_info = Option_Info();
  option_info.bsname = const_cast<char*>(solver_.long_name());
  option_info.wantsol = solver_.wantsol();
  write_solf_ASL(reinterpret_cast<ASL*>(problem_.asl_),
      const_cast<char*>(message.c_str()), const_cast<double*>(values),
      const_cast<double*>(dual_values), &option_info, w.c_str());
}

void mp::SolutionWriter::HandleSolution(
    fmt::StringRef message, const double *values,
    const double *dual_values, double) {
  if (solver_.need_multiple_solutions()) {
    Suffix nsol_suffix = problem_.FindSuffix("nsol", ASL_Sufkind_prob);
    if (nsol_suffix)
      nsol_suffix.set_values(&num_solutions_);
  }
  // TODO: pass to WriteSol
  //option_info.bsname = const_cast<char*>(solver_.long_name());
  //option_info.wantsol = solver_.wantsol();
  const fint *options = problem_.asl_->i.ampl_options_;
  ASLSolutionAdapter sol(problem_,
      message.c_str(), MakeArrayRef(options + 1, options[0]),
      MakeArrayRef(values, values ? problem_.num_vars() : 0),
      MakeArrayRef(dual_values, dual_values ? problem_.num_vars() : 0));
  WriteSol(fmt::format("{}.sol", problem_.name()), sol);
}

mp::ASLSolver::ASLSolver(
    fmt::StringRef name, fmt::StringRef long_name, long date, int flags)
  : SolverImpl<internal::ASLBuilder>(name, long_name, date, flags) {
}

void mp::ASLSolver::RegisterSuffixes(Problem &p) {
  std::size_t num_suffixes = suffixes_.size();
  std::vector<SufDecl> suffix_decls(num_suffixes);
  for (std::size_t i = 0; i < num_suffixes; ++i) {
    const SuffixInfo &si = suffixes_[i];
    SufDecl &sd = suffix_decls[i];
    sd.name = const_cast<char*>(si.name);
    sd.table = const_cast<char*>(si.table);
    sd.kind = si.kind;
    sd.nextra = si.nextra;
  }
  ASL *asl = reinterpret_cast<ASL*>(p.asl_);
  if (asl->i.nsuffixes == 0 && num_suffixes != 0)
    suf_declare_ASL(asl, &suffix_decls[0], static_cast<int>(num_suffixes));
}

bool mp::ASLSolver::ProcessArgs(char **&argv, Problem &p, unsigned flags) {
  ASL *asl = reinterpret_cast<ASL*>(p.asl_);
  RegisterSuffixes(p);

  Option_Info option_info = Option_Info();
  keyword cl_option = keyword();  // command-line option '='
  cl_option.name = const_cast<char*>("=");
  cl_option.desc = const_cast<char*>("show name= possibilities");
  struct OptionPrinter {
    static char *PrintOptionsAndExit(Option_Info *, keyword *kw, char *) {
      Solver *solver = static_cast<Solver*>(kw->info);
      //solver->ShowOptions();
      exit(0);
      return 0;
    }
  };
  cl_option.kf = OptionPrinter::PrintOptionsAndExit;
  cl_option.info = this;
  option_info.options = &cl_option;
  option_info.n_options = 1;

  option_info.sname = const_cast<char*>(name_.c_str());
  option_info.bsname = const_cast<char*>(long_name_.c_str());
  std::string options_var_name = name_ + "_options";
  option_info.opname = const_cast<char*>(options_var_name.c_str());
  option_info.version = const_cast<char*>(version_.c_str());;
  option_info.driver_date = date_;

  char *stub = getstub_ASL(asl, &argv, &option_info);
  if (!stub) {
    PrintUsage(*this, USE_STDERR | FUNC_OPTIONS);
    return false;
  }
  steady_clock::time_point start = steady_clock::now();
  p.Read(stub, read_flags_);
  double read_time = GetTimeAndReset(start);
  bool result = ParseOptions(argv, flags, &p);
  option_info.flags = flags_;
  option_info.wantsol = wantsol_;
  if (timing_)
    Print("Input time = {:.6f}s\n", read_time);
  return result;
}

void mp::ASLSolver::Solve(Problem &p, SolutionHandler &sh) {
  RegisterSuffixes(p);
  DoSolve(p, sh);
}

int mp::ASLSolver::Run(char **argv) {
  Problem p;
  if (!ProcessArgs(argv, p))
    return 1;
  SolutionWriter sol_writer(*this, p);
  Solve(p, sol_writer);
  return 0;
}
