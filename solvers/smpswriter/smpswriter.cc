/*
 SMPS writer implemented as an AMPL solver.

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

#include "smpswriter/smpswriter.h"

#include <cmath>
#include <cstdio>
#include <map>

namespace {

template <typename Map>
typename Map::mapped_type &FindOrInsert(
    Map &map, const typename Map::key_type &key,
    const typename Map::mapped_type &value) {
  typename Map::iterator lb = map.lower_bound(key);
  if (lb != map.end() && !(map.key_comp()(key, lb->first)))
    return lb->second;
  return map.insert(lb, std::make_pair(key, value))->second;
}

std::string ExtractScenario(std::string &name, bool require_scenario = true) {
  std::size_t index_pos = name.find_last_of("[,");
  if (index_pos == std::string::npos) {
    if (!require_scenario)
      return std::string();
    throw mp::Error("Missing scenario index for {}", name);
  }
  bool single_index = name[index_pos] == '[';
  ++index_pos;
  std::string index = name.substr(index_pos, name.size() - index_pos - 1);
  name = name.substr(0, index_pos - 1);
  if (!single_index)
    name += ']';
  return index;
}

double GetConRHSAndType(const mp::Problem &p, int con_index, char &type) {
  mp::Problem::AlgebraicCon con = p.algebraic_con(con_index);
  double lb = con.lb(), ub = con.ub();
  double inf = std::numeric_limits<double>::infinity();
  if (lb <= -inf) {
    type = ub >= inf ? 'N' : 'L';
    return ub;
  }
  if (ub >= inf)
    type = 'G';
  else if (lb == ub)
    type = 'E';
  else
    throw mp::Error("SMPS writer doesn't support ranges");
  return lb;
}
}  // namespace

namespace mp {

class SMPSNameReader {
 private:
  int num_names_;
  internal::NameReader reader_;
  std::vector<fmt::StringRef> names_;

 public:
  explicit SMPSNameReader(int num_names) : num_names_(num_names) {
    names_.reserve(num_names);
  }

  fmt::StringRef name(int index) const { return names_[index]; }

  void OnName(fmt::StringRef name) { names_.push_back(name); }

  void Read(fmt::CStringRef filename) {
    reader_.Read(filename, *this);
    std::size_t num_names = num_names_;
    if (names_.size() != num_names)
      throw mp::Error("expected {} names in {}", num_names_, filename);
  }
};

class FileWriter {
 private:
  FILE *f_;
  FMT_DISALLOW_COPY_AND_ASSIGN(FileWriter);

 public:
  explicit FileWriter(fmt::CStringRef name)
    : f_(std::fopen(name.c_str(), "w")) {}
  ~FileWriter() { std::fclose(f_); }

  void Write(fmt::CStringRef format, const fmt::ArgList &args) {
    fmt::print(f_, format, args);
  }
  FMT_VARIADIC(void, Write, fmt::CStringRef)
};

SMPSWriter::SMPSWriter()
  : SolverImpl<ColProblem>("smpswriter", "SMPSWriter", 20130709) {
  AddSuffix("stage", 0, suf::VAR);
}

void SMPSWriter::SplitConRHSIntoScenarios(
    const Problem &p, std::vector<CoreConInfo> &core_cons,
    SMPSNameReader &con_names) {
  int num_cons = p.num_algebraic_cons();
  for (int i = 0; i < num_cons; ++i) {
    CoreConInfo &info = core_cons[con_info[i].core_index];
    int scenario_index = con_info[i].scenario_index;
    if (scenario_index != 0 && info.type)
      continue;
    double rhs = GetConRHSAndType(p, i, info.type);
    if (scenario_index == 0)
      info.rhs = rhs;
  }
  for (int i = 0; i < num_cons; ++i) {
    int scenario_index = con_info[i].scenario_index;
    if (scenario_index == 0)
      continue;
    char type = 0;
    double rhs = GetConRHSAndType(p, i, type);
    int core_con_index = con_info[i].core_index;
    const CoreConInfo &info = core_cons[core_con_index];
    if (type != info.type) {
      throw Error("Inconsistent constraint type for {}", con_names.name(i));
    }
    if (rhs != info.rhs)
      scenarios[scenario_index].AddRHS(core_con_index, rhs);
  }
}

void SMPSWriter::SplitVarBoundsIntoScenarios(
    const Problem &p, std::vector<CoreVarInfo> &core_vars) {
  int num_vars = p.num_vars();
  for (int i = 0; i < num_vars; ++i) {
    int scenario_index = var_info[i].scenario_index;
    if (scenario_index != 0)
      continue;
    CoreVarInfo &info = core_vars[var_info[i].core_index];
    Problem::Variable var = p.var(i);
    info.lb = var.lb();
    info.ub = var.ub();
  }
  for (int i = 0; i < num_vars; ++i) {
    int scenario_index = var_info[i].scenario_index;
    if (scenario_index == 0)
      continue;
    int core_var_index = var_info[i].core_index;
    const CoreVarInfo &info = core_vars[core_var_index];
    Problem::Variable var = p.var(i);
    double lb = var.lb(), ub = var.ub();
    if (lb != info.lb)
      scenarios[scenario_index].AddLB(core_var_index, lb);
    if (ub != info.ub)
      scenarios[scenario_index].AddUB(core_var_index, ub);
  }
}

void SMPSWriter::WriteColumns(
    FileWriter &writer, const ColProblem &p, int num_stages,
    int num_core_cons, const std::vector<double> &core_obj_coefs) {
  writer.Write("COLUMNS\n");
  std::vector<double> core_coefs(num_core_cons);
  std::vector<int> nonzero_coef_indices;
  nonzero_coef_indices.reserve(num_core_cons);
  int int_var_index = 0;
  IntSuffix stage_suffix = p.suffixes(suf::VAR).Find<int>("stage");
  bool integer_block = false;
  for (int stage = 0; stage < num_stages; ++stage) {
    for (int i = 0, n = p.num_vars(); i < n; ++i) {
      int var_stage = stage_suffix ? std::max(stage_suffix.value(i) - 1, 0) : 0;
      if (var_stage != stage) continue;
      int core_var_index = var_info[i].core_index;
      if (var_info[i].scenario_index == 0) {
        // Clear the core_coefs vector.
        for (std::vector<int>::const_iterator j = nonzero_coef_indices.begin(),
            end = nonzero_coef_indices.end(); j != end; ++j) {
          core_coefs[*j] = 0;
        }
        nonzero_coef_indices.clear();

        if (p.var(i).type() == var::CONTINUOUS) {
          if (integer_block) {
            writer.Write(
                "    INT{:<5}    'MARKER'      'INTEND'\n", int_var_index);
            integer_block = false;
          }
        } else if (!integer_block) {
          writer.Write(
              "    INT{:<5}    'MARKER'      'INTORG'\n", ++int_var_index);
          integer_block = true;
        }

        if (double obj_coef = core_obj_coefs[core_var_index]) {
          writer.Write(
              "    C{:<7}  OBJ       {}\n", core_var_index + 1, obj_coef);
        }

        // Write the core coefficients and store them in the core_coefs vector.
        for (int k = p.col_start(i), end = p.col_start(i + 1); k != end; ++k) {
          int con_index = p.row_index(k);
          if (con_info[con_index].scenario_index != 0)
            continue;
          int core_con_index = con_info[con_index].core_index;
          core_coefs[core_con_index] = p.value(k);
          nonzero_coef_indices.push_back(core_con_index);
          writer.Write("    C{:<7}  R{:<7}  {}\n",
              core_var_index + 1, core_con_index + 1, p.value(k));
        }
      }

      // Go over non-core coefficients and compare them to those in the core.
      for (int k = p.col_start(i), end = p.col_start(i + 1); k != end; ++k) {
        int con_index = p.row_index(k);
        int scenario_index = con_info[con_index].scenario_index;
        if (scenario_index == 0)
          continue;
        int core_con_index = con_info[con_index].core_index;
        double coef = p.value(k);
        double core_coef = core_coefs[core_con_index];
        if (coef != core_coef) {
          scenarios[scenario_index].AddConTerm(
              core_con_index, core_var_index, coef);
          if (core_coef == 0) {
            writer.Write("    C{:<7}  R{:<7}  0\n",
                core_var_index + 1, core_con_index + 1);
          }
        }
      }
    }
  }
  if (integer_block)
    writer.Write("    INT{:<5}    'MARKER'      'INTEND'\n", int_var_index);
}

void SMPSWriter::Solve(ColProblem &p, SolutionHandler &) {
  if (p.has_nonlinear_cons())
    throw Error("SMPS writer doesn't support nonlinear problems");

  // Count the number of stages and the number of variables in stage 0.
  int num_vars = p.num_vars();
  int num_stage0_vars = num_vars;
  int num_stages = 1;
  IntSuffix stage_suffix = p.suffixes(suf::VAR).Find<int>("stage");
  if (stage_suffix) {
    num_stage0_vars = 0;
    for (int i = 0; i < num_vars; ++i) {
      int stage_plus_1 = stage_suffix.value(i);
      if (stage_plus_1 > 0)
        num_stages = std::max(stage_plus_1, num_stages);
      else
        ++num_stage0_vars;
    }
  }
  if (num_stages > 2)
    throw Error("SMPS writer doesn't support problems with more than 2 stages");

  std::string smps_basename = p.name();
  std::string::size_type ext_pos = smps_basename.rfind('.');
  if (ext_pos != std::string::npos)
    smps_basename.resize(ext_pos);

  int num_cons = p.num_algebraic_cons();
  int num_stage0_cons = num_cons;
  int num_core_vars = num_vars, num_core_cons = num_cons;
  var_info.resize(num_vars);
  con_info.resize(num_cons);
  scenarios.clear();
  SMPSNameReader var_names(p.num_vars());
  SMPSNameReader con_names(p.num_algebraic_cons() + p.num_objs());
  if (stage_suffix) {
    var_names.Read(smps_basename + ".col");
    con_names.Read(smps_basename + ".row");
    std::map<std::string, int> scenario_indices;
    int stage0_var_count = 0;
    std::map<std::string, int> stage1_vars;
    for (int i = 0; i < num_vars; ++i) {
      int stage = stage_suffix.value(i) - 1;
      VarConInfo &info = var_info[i];
      if (stage > 0) {
        // Split the name into scenario and the rest and merge variables that
        // only differ by scenario into the same variable.
        std::string name = var_names.name(i).to_string();
        std::string scenario = ExtractScenario(name);
        info.core_index = FindOrInsert(stage1_vars, name,
          static_cast<int>(num_stage0_vars + stage1_vars.size()));
        info.scenario_index = FindOrInsert(scenario_indices,
          scenario, static_cast<int>(scenario_indices.size()));
      } else {
        info.core_index = stage0_var_count++;
      }
    }
    assert(stage0_var_count == num_stage0_vars);
    num_core_vars = static_cast<int>(num_stage0_vars + stage1_vars.size());

    // Compute stage of each constraint as a maximum of stages of
    // variables in it.
    std::vector<int> con_stages(num_cons);
    std::map<std::string, int> stage1_cons;
    for (int j = 0; j < num_vars; ++j) {
      int stage = stage_suffix.value(j) - 1;
      if (stage <= 0) continue;
      // Update stages of all constraints containing this variable.
      for (int k = p.col_start(j), end = p.col_start(j + 1); k != end; ++k) {
        int con_index = p.row_index(k);
        int &con_stage = con_stages[con_index];
        if (stage > con_stage) {
          if (con_stage == 0) {
            // Split the name into scenario and the rest and merge constraints
            // that only differ by scenario into the same constraint.
            std::string name = con_names.name(con_index).to_string();
            std::string scenario = ExtractScenario(name);
            VarConInfo &info = con_info[con_index];
            info.core_index = FindOrInsert(
                stage1_cons, name, static_cast<int>(stage1_cons.size()));
            info.scenario_index = FindOrInsert(scenario_indices,
                scenario, static_cast<int>(scenario_indices.size()));
            --num_stage0_cons;
          }
          con_stage = stage;
        }
      }
    }

    // Check for any stage 1 constraints remaining - some of them may not
    // have been detected previously because they don't have stage 1 variables.
    // This can happen if all stage 1 variables have zero coefficients in a
    // core constraint.
    int stage0_con_count = 0;
    for (int i = 0; i < num_cons; ++i) {
      if (con_stages[i] != 0)
        continue;
      VarConInfo &info = con_info[i];
      // A constraint with a name which only differs in scenario from a
      // name of some other stage 1 constraint, is also a stage 1 constraint.
      std::string name = con_names.name(i).to_string();
      std::string scenario = ExtractScenario(name, false);
      std::map<std::string, int>::iterator stage1_con =
          scenario.empty() ? stage1_cons.end() : stage1_cons.find(name);
      if (stage1_con != stage1_cons.end()) {
        con_stages[i] = 1;
        info.core_index = stage1_con->second;
        info.scenario_index = FindOrInsert(scenario_indices,
            scenario, static_cast<int>(scenario_indices.size()));
        --num_stage0_cons;
      } else {
        info.core_index = stage0_con_count++;
      }
    }
    assert(stage0_con_count == num_stage0_cons);
    num_core_cons = static_cast<int>(num_stage0_cons + stage1_cons.size());

    for (int i = 0; i < num_cons; ++i) {
      if (con_stages[i] != 0)
        con_info[i].core_index += num_stage0_cons;
    }

    scenarios.resize(scenario_indices.size());
  } else {
    for (int i = 0; i < num_vars; ++i)
      var_info[i].core_index = i;
    for (int i = 0; i < num_cons; ++i)
      con_info[i].core_index = i;
    scenarios.resize(1);
  }

  std::vector<CoreConInfo> core_cons(num_core_cons);
  SplitConRHSIntoScenarios(p, core_cons, con_names);
  std::vector<CoreVarInfo> core_vars(num_core_vars);
  SplitVarBoundsIntoScenarios(p, core_vars);

  // Write the .tim file.
  {
    FileWriter writer(smps_basename + ".tim");
    writer.Write(
      "TIME          PROBLEM\n"
      "PERIODS\n"
      "    C1        OBJ                      T1\n");
    if (num_stages > 1) {
      writer.Write("    C{:<7}  R{:<7}                 T2\n",
          num_stage0_vars + 1, num_stage0_cons + 1);
    }
    writer.Write("ENDATA\n");
  }

  // Write the .cor file.
  std::vector<double> probabilities(scenarios.size());
  {
    FileWriter writer(smps_basename + ".cor");
    writer.Write(
      "NAME          PROBLEM\n"
      "ROWS\n"
      " N  OBJ\n");
    for (int i = 0; i < num_core_cons; ++i)
      writer.Write(" {}  R{}\n", core_cons[i].type, i + 1);

    std::vector<double> core_obj_coefs(num_core_vars);
    std::vector<double> sum_core_obj_coefs;
    if (p.num_objs() != 0) {
      Problem::Objective obj = p.obj(0);
      const LinearExpr &obj_expr = obj.linear_expr();
      if (obj.nonlinear_expr())
        throw Error("SMPS writer doesn't support nonlinear problems");
      int reference_var_index = 0;
      int core_reference_var_index = 0;
      if (probabilities.size() != 1) {
        // Deduce probabilities from objective coefficients.
        for (LinearExpr::iterator
             i = obj_expr.begin(), end = obj_expr.end(); i != end; ++i) {
          int stage = stage_suffix.value(i->var_index()) - 1;
          if (stage > 0) {
            reference_var_index = i->var_index();
            core_reference_var_index = var_info[i->var_index()].core_index;
            break;
          }
        }
        sum_core_obj_coefs.resize(num_core_vars);
        for (LinearExpr::iterator
             i = obj_expr.begin(), end = obj_expr.end(); i != end; ++i) {
          const VarConInfo &info = var_info[i->var_index()];
          if (info.core_index == core_reference_var_index)
            probabilities[info.scenario_index] = i->coef();
          sum_core_obj_coefs[info.core_index] += i->coef();
        }
        for (size_t i = 0, n = scenarios.size(); i != n; ++i)
          probabilities[i] /= sum_core_obj_coefs[core_reference_var_index];
      } else {
        probabilities[0] = 1;
      }

      // Compute objective coefficients in the core problem.
      for (LinearExpr::iterator
           i = obj_expr.begin(), end = obj_expr.end(); i != end; ++i) {
        const VarConInfo &info = var_info[i->var_index()];
        if (info.scenario_index == 0) {
          double coef = i->coef();
          if (stage_suffix && stage_suffix.value(i->var_index()) - 1 > 0)
            coef /= probabilities[0];
          core_obj_coefs[info.core_index] = coef;
        }
        // Check probabilities deduced using other variables.
        if (probabilities.size() != 1 && stage_suffix &&
            stage_suffix.value(i->var_index()) - 1 > 0) {
          double ref_prob = probabilities[info.scenario_index];
          double prob = i->coef() / sum_core_obj_coefs[info.core_index];
          double prob_tolerance = 1e-5;
          if (std::abs(prob - ref_prob) > prob_tolerance) {
            throw Error("Probability deduced using variable {} ({}) "
                "is inconsistent with the one deduced using variable {} ({})",
                    var_names.name(reference_var_index),
                    probabilities[info.scenario_index],
                    var_names.name(i->var_index()), prob);
          }
        }
      }
    }

    WriteColumns(writer, p, num_stages, num_core_cons, core_obj_coefs);

    writer.Write("RHS\n");
    for (int i = 0; i < num_core_cons; ++i)
      writer.Write("    RHS1      R{:<7}  {}\n", i + 1, core_cons[i].rhs);

    writer.Write("BOUNDS\n");
    double inf = std::numeric_limits<double>::infinity();
    for (int i = 0; i < num_core_vars; ++i) {
      double lb = core_vars[i].lb, ub = core_vars[i].ub;
      if (lb != 0)
        writer.Write(" LO BOUND1      C{:<7}  {}\n", i + 1, lb);
      if (ub < inf)
        writer.Write(" UP BOUND1      C{:<7}  {}\n", i + 1, ub);
    }

    writer.Write("ENDATA\n");
  }

  // Write the .sto file.
  {
    FileWriter writer(smps_basename + ".sto");
    writer.Write(
      "STOCH         PROBLEM\n"
      "SCENARIOS     DISCRETE\n");
    if (num_stages > 1) {
      writer.Write(" SC SCEN1     'ROOT'    {:<12}   T1\n", probabilities[0]);
      for (size_t i = 1, n = scenarios.size(); i < n; ++i) {
        writer.Write(" SC SCEN{:<4}  SCEN1     {:<12}   T2\n",
            i + 1, probabilities[i]);
        for (Scenario::ConTermIterator j = scenarios[i].con_term_begin(),
            end = scenarios[i].con_term_end(); j != end; ++j) {
          writer.Write("    C{:<7}  R{:<7}  {}\n",
              j->var_index + 1, j->con_index + 1, j->coef);
        }
        for (Scenario::RHSIterator j = scenarios[i].rhs_begin(),
            end = scenarios[i].rhs_end(); j != end; ++j) {
          writer.Write("    RHS1      R{:<7}  {}\n",
              j->con_index + 1, j->rhs);
        }
        for (Scenario::BoundIterator j = scenarios[i].lb_begin(),
            end = scenarios[i].lb_end(); j != end; ++j) {
          writer.Write(" LO BOUND1      C{:<7}  {}\n",
              j->var_index + 1, j->bound);
        }
        for (Scenario::BoundIterator j = scenarios[i].ub_begin(),
            end = scenarios[i].ub_end(); j != end; ++j) {
          writer.Write(" UP BOUND1      C{:<7}  {}\n",
              j->var_index + 1, j->bound);
        }
      }
    }
    writer.Write("ENDATA\n");
  }
}

SolverPtr create_smpswriter(const char *) {
  return SolverPtr(new SMPSWriter());
}
}  // namespace mp
