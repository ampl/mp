/*
 IBM/ILOG CPLEX solver for AMPL.

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

#include "solvers/ilogcp/cplex.h"

#include <cctype>
#include <cstdlib>
#include <vector>

#include "solvers/util/clock.h"
#include "solvers/ilogcp/concert.h"
#include "solvers/ilogcp/ilogcp_date.h"

using std::vector;

namespace {

ampl::OptionError GetOptionValueError(
    fmt::StringRef name, fmt::StringRef message) {
  throw ampl::OptionError(fmt::Format(
      "Can't get value of option {}: {}") << name.c_str() << message.c_str());
}
}

namespace ampl {

std::string CPLEXSolver::GetOptionHeader() {
  return "IloCPLEX Directives for AMPL\n"
      "--------------------------\n"
      "\n"
      "To set these directives, assign a string specifying their values to the AMPL "
      "option ilocplex_options.  For example:\n"
      "\n"
      "  ampl: option ilocplex_options 'mipdisplay=2 mipinterval=10';\n"
      "\n"
      "Where both a number and a keyword are given, either may be used to specify "
      "the option setting.\n";
}

CPLEXSolver::CPLEXSolver() :
   Solver("ilocplex", 0, YYYYMMDD), cplex_(env_), aborter_(env_) {
  options_[DEBUGEXPR] = false;
  options_[USENUMBEROF] = true;

  set_long_name(fmt::Format("ilocplex {}.{}.{}")
      << IloConcertVersion::_ILO_MAJOR_VERSION
      << IloConcertVersion::_ILO_MINOR_VERSION
      << IloConcertVersion::_ILO_TECH_VERSION);
  set_version(fmt::Format("AMPL/IBM ILOG CPLEX Optimizer [{} {}.{}.{}]")
      << IloConcertVersion::_ILO_NAME << IloConcertVersion::_ILO_MAJOR_VERSION
      << IloConcertVersion::_ILO_MINOR_VERSION
      << IloConcertVersion::_ILO_TECH_VERSION);

  AddIntOption("debugexpr",
      "0 or 1 (default 0):  Whether to print debugging "
      "information for expression trees.",
      &CPLEXSolver::GetBoolOption, &CPLEXSolver::SetBoolOption, DEBUGEXPR);

  AddIntOption<CPLEXSolver, int>("mipdisplay",
      "Frequency of displaying branch-and-bound information "
      "(for optimizing integer variables):\n"
      "      0 (default) = never\n"
      "      1 = each integer feasible solution\n"
      "      2 = every \"mipinterval\" nodes\n"
      "      3 = every \"mipinterval\" nodes plus\n"
      "          information on LP relaxations\n"
      "          (as controlled by \"display\")\n"
      "      4 = same as 2, plus LP relaxation info.\n"
      "      5 = same as 2, plus LP subproblem info.\n",
      &CPLEXSolver::GetCPLEXIntOption, &CPLEXSolver::SetCPLEXIntOption,
      IloCplex::MIPDisplay);

  AddIntOption<CPLEXSolver, int>("mipinterval",
      "Frequency of node logging for mipdisplay 2 or 3. Default = 1.",
      &CPLEXSolver::GetCPLEXIntOption, &CPLEXSolver::SetCPLEXIntOption,
      IloCplex::MIPInterval);

  AddIntOption("usenumberof",
      "0 or 1 (default 1):  Whether to consolidate 'numberof' expressions "
      "by use of IloDistribute constraints.",
      &CPLEXSolver::GetBoolOption, &CPLEXSolver::SetBoolOption,
      CPLEXSolver::USENUMBEROF);
}

CPLEXSolver::~CPLEXSolver() {
  env_.end();
}

void CPLEXSolver::SetBoolOption(const char *name, int value, Option opt) {
  if (value != 0 && value != 1)
    throw InvalidOptionValue(name, value);
  options_[opt] = value;
}

int CPLEXSolver::GetCPLEXIntOption(const char *name, int param) const {
  // Use CPXgetintparam instead of IloCplex::setParam to avoid dealing with
  // two overloads, one for the type int and one for the type long.
  int value = 0;
  int result = CPXgetintparam(cplex_.getImpl()->getCplexEnv(), param, &value);
  if (result != 0)
    throw GetOptionValueError(name, fmt::Format("CPLEX error = {}") << result);
  return value;
}

void CPLEXSolver::SetCPLEXIntOption(const char *name, int value, int param) {
  // Use CPXsetintparam instead of IloCplex::setParam to avoid dealing with
  // two overloads, one for the type int and one for the type long.
  if (CPXsetintparam(cplex_.getImpl()->getCplexEnv(), param, value) != 0)
    throw InvalidOptionValue(name, value);
}

void CPLEXSolver::DoSolve(Problem &p) {
  steady_clock::time_point time = steady_clock::now();

  NLToConcertConverter converter(env_,
      GetOption(USENUMBEROF), GetOption(DEBUGEXPR));
  converter.Convert(p);
  IloModel model = converter.model();
  IloNumVarArray vars = converter.vars();

  try {
    cplex_.extract(model);
  } catch (IloAlgorithm::CannotExtractException &e) {
    const IloExtractableArray &extractables = e.getExtractables();
    if (extractables.getSize() == 0)
      throw;
    throw UnsupportedExprError::CreateFromExprString(
        str(fmt::Format("{}") << extractables[0]));
  }
  cplex_.setParam(IloCplex::MIPDisplay, 0);
  cplex_.use(aborter_);
  SignalHandler sh(*this, this);

  double setup_time = GetTimeAndReset(time);
  cplex_.solve();
  double solution_time = GetTimeAndReset(time);

  // Convert solution status.
  int solve_code = 0;
  bool has_solution = false;
  std::string status =
      ConvertSolutionStatus(cplex_, sh, solve_code, has_solution);
  p.set_solve_code(solve_code);

  fmt::Writer writer;
  writer.Format("{}: {}\n") << long_name() << status;
  double obj_value = std::numeric_limits<double>::quiet_NaN();
  vector<double> solution, dual_solution;
  if (has_solution) {
    int num_vars = p.num_vars();
    solution.resize(num_vars);
    for (int j = 0; j < num_vars; ++j) {
      IloNumVar &v = vars[j];
      solution[j] = cplex_.isExtracted(v) ? cplex_.getValue(v) : v.getLB();
    }

    if (cplex_.isMIP()) {
      writer << cplex_.getNnodes() << " nodes, ";
    } else {
      IloRangeArray cons = converter.cons();
      IloInt num_cons = cons.getSize();
      dual_solution.resize(num_cons);
      for (IloInt i = 0; i < num_cons; ++i)
        dual_solution[i] = cplex_.getDual(cons[i]);
    }
    writer << cplex_.getNiterations() << " iterations";

    if (p.num_objs() > 0) {
      obj_value = cplex_.getObjValue();
      writer.Format(", objective {}") << ObjPrec(obj_value);
    }
  }
  HandleSolution(p, writer.c_str(),
      ptr(solution), ptr(dual_solution), obj_value);
  double output_time = GetTimeAndReset(time);

  if (timing()) {
    Print("Setup time = {:.6f}s\n"
          "Solution time = {:.6f}s\n"
          "Output time = {:.6f}s\n")
            << setup_time << solution_time << output_time;
  }
}
}
