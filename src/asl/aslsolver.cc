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
#undef Char
}

#include "mp/clock.h"
#include "mp/solver.h"
#include "problem.h"

void mp::Solver::RegisterSuffixes(Problem &p) {
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

void mp::Solver::SolutionWriter::HandleFeasibleSolution(
    Problem &p, fmt::StringRef message, const double *values,
    const double *dual_values, double) {
  ++solver_->num_solutions_;
  if (solver_->solution_stub_.empty())
    return;
  fmt::Writer w;
  w << solver_->solution_stub_ << solver_->num_solutions_ << ".sol";
  Option_Info option_info = Option_Info();
  option_info.bsname = const_cast<char*>(solver_->long_name_.c_str());
  option_info.wantsol = solver_->wantsol_;
  write_solf_ASL(reinterpret_cast<ASL*>(p.asl_),
      const_cast<char*>(message.c_str()), const_cast<double*>(values),
      const_cast<double*>(dual_values), &option_info, w.c_str());
}

void mp::Solver::SolutionWriter::HandleSolution(
    Problem &p, fmt::StringRef message, const double *values,
    const double *dual_values, double) {
  if (solver_->need_multiple_solutions()) {
    Suffix nsol_suffix = p.suffix("nsol", ASL_Sufkind_prob);
    if (nsol_suffix)
      nsol_suffix.set_values(&solver_->num_solutions_);
  }
  Option_Info option_info = Option_Info();
  option_info.bsname = const_cast<char*>(solver_->long_name_.c_str());
  option_info.wantsol = solver_->wantsol_;
  write_sol_ASL(reinterpret_cast<ASL*>(p.asl_),
      const_cast<char*>(message.c_str()), const_cast<double*>(values),
      const_cast<double*>(dual_values), &option_info);
}

bool mp::Solver::ProcessArgs(char **&argv, Problem &p, unsigned flags) {
  ASL *asl = reinterpret_cast<ASL*>(p.asl_);
  RegisterSuffixes(p);

  // Workaround for GCC bug 30111 that prevents value-initialization of
  // the base POD class.
  Option_Info option_info = Option_Info();
  keyword cl_option = keyword();  // command-line option '='
  cl_option.name = const_cast<char*>("=");
  cl_option.desc = const_cast<char*>("show name= possibilities");
  struct OptionPrinter {
    static char *PrintOptionsAndExit(Option_Info *, keyword *kw, char *) {
      Solver *solver = static_cast<Solver*>(kw->info);
      fmt::Writer writer;
      internal::FormatRST(writer, solver->option_header_);
      if (!solver->option_header_.empty())
        writer << '\n';
      writer << "Options:\n";
      const int DESC_INDENT = 6;
      const OptionSet &options = solver->options_;
      for (OptionSet::const_iterator
           i = options.begin(); i != options.end(); ++i) {
        SolverOption *opt = *i;
        writer << '\n' << opt->name() << '\n';
        internal::FormatRST(writer, opt->description(),
                            DESC_INDENT, opt->values());
      }
      std::fwrite(writer.data(), writer.size(), 1, stdout);
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

  // Make "-?" show some command-line options that are relevant when
  // importing functions.
  option_info.flags |= ASL_OI_want_funcadd;

  char *stub = getstub_ASL(asl, &argv, &option_info);
  if (!stub) {
    usage_noexit_ASL(&option_info, 1);
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

int mp::Solver::Run(char **argv) {
  Problem p;
  if (!ProcessArgs(argv, p))
    return 1;
  Solve(p);
  return 0;
}
