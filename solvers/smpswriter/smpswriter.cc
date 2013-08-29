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

#include "smpswriter.h"

#include <cstdio>
#include <set>

namespace {
SufDecl STAGE_SUFFIX = {const_cast<char*>("stage"), 0, ASL_Sufkind_var};

std::string ExtractScenario(std::string &name, bool require_scenario = true) {
  std::size_t index_pos = name.find_last_of("[,");
  if (index_pos == std::string::npos) {
    if (!require_scenario)
      return std::string();
    throw ampl::Error(fmt::Format("Missing scenario index for {}") << name);
  }
  bool single_index = name[index_pos] == '[';
  ++index_pos;
  std::string index = name.substr(index_pos, name.size() - index_pos - 1);
  name = name.substr(0, index_pos - 1);
  if (!single_index)
    name += ']';
  return index;
}
}

namespace ampl {

SMPSWriter::SMPSWriter() :
  Solver<SMPSWriter>("smpswriter", "SMPSWriter", 20130709) {
  DeclareSuffixes(&STAGE_SUFFIX, 1);
  set_read_flags(Problem::READ_COLUMNWISE);
}

class FileWriter : Noncopyable {
 private:
  FILE *f_;

  class Writer {
   private:
    FILE *f_;

   public:
    Writer(FILE *f) : f_(f) {}

    void operator()(const fmt::BasicFormatter<char> &f) const {
      std::fwrite(f.data(), 1, f.size(), f_);
    }
  };

 public:
  FileWriter(fmt::StringRef name) : f_(std::fopen(name.c_str(), "w")) {}
  ~FileWriter() { std::fclose(f_); }

  fmt::TempFormatter<Writer> Write(fmt::StringRef format) {
    return fmt::TempFormatter<Writer>(format, Writer(f_));
  }
};

template <typename Map>
typename Map::mapped_type &FindOrInsert(
    Map &map, const typename Map::key_type &key,
    const typename Map::mapped_type &value) {
  auto lb = map.lower_bound(key);
  if (lb != map.end() && !(map.key_comp()(key, lb->first)))
    return lb->second;
  return map.insert(lb, std::make_pair(key, value))->second;
}

class Scenario {
 public:
  struct ConTerm {
    int con_index;
    int var_index;
    double coef;
    ConTerm(int con_index, int var_index, double coef)
    : con_index(con_index), var_index(var_index), coef(coef) { }
  };

 private:
  std::vector<ConTerm> con_terms_;

 public:
  Scenario() {}

  void Add(int con_index, int var_index, double coef) {
    con_terms_.push_back(ConTerm(con_index, var_index, coef));
  }

  typedef std::vector<ConTerm>::const_iterator iterator;

  iterator begin() const { return con_terms_.begin(); }
  iterator end() const { return con_terms_.end(); }
};

void SMPSWriter::Solve(Problem &p) {
  if (p.num_nonlinear_objs() != 0 || p.num_nonlinear_cons() != 0)
    throw Error("SMPS writer doesn't support nonlinear problems");

  // Count the number of stages and the number of variables in stage 0.
  int num_vars = p.num_vars();
  Suffix stage_suffix = p.suffix("stage", ASL_Sufkind_var);
  int num_stage0_vars = num_vars;
  int num_stages = 1;
  if (stage_suffix) {
    num_stage0_vars = 0;
    for (int i = 0; i < num_vars; ++i) {
      int stage_plus_1 = stage_suffix.int_value(i);
      if (stage_plus_1 > 0)
        num_stages = std::max(stage_plus_1, num_stages);
      else
        ++num_stage0_vars;
    }
  }
  if (num_stages > 2)
    throw Error("SMPS writer doesn't support problems with more than 2 stages");

  // Information about a variable or constraint.
  struct VarConInfo {
    int core_index;  // index of this variable in the core problem
    int scenario_index;
    VarConInfo() : core_index(0), scenario_index() {}
  };

  int num_cons = p.num_cons();
  int num_stage0_cons = num_cons;
  int num_core_vars = num_vars, num_core_cons = num_cons;
  std::vector<VarConInfo> var_info(num_vars);
  std::vector<VarConInfo> con_info(num_cons);
  std::vector<Scenario> scenarios;
  if (stage_suffix) {
    std::map<std::string, int> scenario_indices;
    int stage0_var_count = 0;
    std::map<std::string, int> stage1_vars;
    for (int i = 0; i < num_vars; ++i) {
      int stage = stage_suffix.int_value(i) - 1;
      auto &info = var_info[i];
      if (stage > 0) {
        // Split the name into scenario and the rest and merge variables that
        // only differ by scenario into the same variable.
        std::string name = p.var_name(i);
        std::string scenario = ExtractScenario(name);
        info.core_index = FindOrInsert(
            stage1_vars, name, num_stage0_vars + stage1_vars.size());
        info.scenario_index =
            FindOrInsert(scenario_indices, scenario, scenario_indices.size());
      } else {
        info.core_index = stage0_var_count++;
      }
    }
    assert(stage0_var_count == num_stage0_vars);
    num_core_vars = num_stage0_vars + stage1_vars.size();

    // Compute stage of each constraint as a maximum of stages of
    // variables in it.
    Problem::ColMatrix matrix = p.col_matrix();
    std::vector<int> con_stages(num_cons);
    std::map<std::string, int> stage1_cons;
    for (int j = 0; j < num_vars; ++j) {
      int stage = stage_suffix.int_value(j) - 1;
      if (stage <= 0) continue;
      // Update stages of all constraints containing this variable.
      for (int k = matrix.col_start(j),
          end = matrix.col_start(j + 1); k != end; ++k) {
        int con_index = matrix.row_index(k);
        int &con_stage = con_stages[con_index];
        if (stage > con_stage) {
          if (con_stage == 0) {
            // Split the name into scenario and the rest and merge constraints
            // that only differ by scenario into the same constraint.
            std::string name = p.con_name(con_index);
            std::string scenario = ExtractScenario(name);
            auto &info = con_info[con_index];
            info.core_index = FindOrInsert(
                stage1_cons, name, stage1_cons.size());
            info.scenario_index = FindOrInsert(
                scenario_indices, scenario, scenario_indices.size());
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
      auto &info = con_info[i];
      // A constraint with a name which only differs in scenario from a
      // name of some other stage 1 constraint, is also a stage 1 constraint.
      std::string name = p.con_name(i);
      std::string scenario = ExtractScenario(name, false);
      auto stage1_con = scenario.empty() ?
          stage1_cons.end() : stage1_cons.find(name);
      if (stage1_con != stage1_cons.end()) {
        con_stages[i] = 1;
        info.core_index = stage1_con->second;
        info.scenario_index = FindOrInsert(
            scenario_indices, scenario, scenario_indices.size());
        --num_stage0_cons;
      } else {
        info.core_index = stage0_con_count++;
      }
    }
    assert(stage0_con_count == num_stage0_cons);
    num_core_cons = num_stage0_cons + stage1_cons.size();

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

  std::string smps_basename = p.name();

  // Write the .tim file.
  {
    FileWriter writer(smps_basename + ".tim");
    writer.Write(
      "TIME          PROBLEM\n"
      "PERIODS\n"
      "    C1        OBJ                      T1\n");
    if (num_stages > 1) {
      writer.Write("    C{:<7}  R{:<7}                 T2\n")
          << num_stage0_vars + 1 << num_stage0_cons + 1;
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
    for (int i = 0; i < num_cons; ++i) {
      if (con_info[i].scenario_index != 0)
        continue;
      double lb = p.con_lb(i), ub = p.con_ub(i);
      char type = 0;
      if (lb <= negInfinity)
        type = ub >= Infinity ? 'N' : 'L';
      else
        type = ub >= Infinity ? 'G' : 'E';
      writer.Write(" {}  R{}\n") << type << con_info[i].core_index + 1;
    }

    writer.Write("COLUMNS\n");
    std::vector<double> core_obj_coefs(num_core_vars);
    if (p.num_objs() != 0) {
      LinearObjExpr obj_expr = p.linear_obj_expr(0);
      if (probabilities.size() != 1) {
        // Deduce probabilities from objective coefficients.
        int reference_var_index = 0;
        for (auto i = obj_expr.begin(), end = obj_expr.end(); i != end; ++i) {
          int stage = stage_suffix.int_value(i->var_index()) - 1;
          if (stage > 0) {
            reference_var_index = var_info[i->var_index()].core_index;
            break;
          }
        }
        std::vector<double> sum_core_obj_coefs(num_core_vars);
        for (auto i = obj_expr.begin(), end = obj_expr.end(); i != end; ++i) {
          const VarConInfo &info = var_info[i->var_index()];
          if (info.core_index == reference_var_index)
            probabilities[info.scenario_index] = i->coef();
          sum_core_obj_coefs[info.core_index] += i->coef();
        }
        for (size_t i = 0, n = scenarios.size(); i != n; ++i)
          probabilities[i] /= sum_core_obj_coefs[reference_var_index];
      } else {
        probabilities[0] = 1;
      }

      // Compute objective coefficients in the core problem.
      for (auto i = obj_expr.begin(), end = obj_expr.end(); i != end; ++i) {
        const VarConInfo &info = var_info[i->var_index()];
        if (info.scenario_index == 0) {
          double coef = i->coef();
          if (stage_suffix && stage_suffix.int_value(i->var_index()) - 1 > 0)
            coef /= probabilities[0];
          core_obj_coefs[info.core_index] = coef;
        }
        // TODO: check probabilities deduced from other variables
      }
    }
    std::vector<double> core_coefs(num_core_cons);
    std::vector<int> nonzero_coef_indices;
    nonzero_coef_indices.reserve(num_core_cons);
    for (int i = 0; i < num_vars; ++i) {
      int core_var_index = var_info[i].core_index;
      Problem::ColMatrix matrix = p.col_matrix();
      // TODO: what if the variables are not ordered by scenarios?
      if (var_info[i].scenario_index == 0) {
        // Clear the core_coefs vector.
        for (auto j = nonzero_coef_indices.begin(),
            end = nonzero_coef_indices.end(); j != end; ++j) {
          core_coefs[*j] = 0;
        }
        nonzero_coef_indices.clear();

        if (double obj_coef = core_obj_coefs[core_var_index]) {
          writer.Write("    C{:<7}  OBJ       {}\n")
              << core_var_index + 1 << obj_coef;
        }

        // Write the core coefficients and store them in the core_coefs vector.
        for (int k = matrix.col_start(i),
            end = matrix.col_start(i + 1); k != end; ++k) {
          int con_index = matrix.row_index(k);
          if (con_info[con_index].scenario_index != 0)
            continue;
          int core_con_index = con_info[con_index].core_index;
          core_coefs[core_con_index] = matrix.value(k);
          nonzero_coef_indices.push_back(core_con_index);
          writer.Write("    C{:<7}  R{:<7}  {}\n")
              << core_var_index + 1 << core_con_index + 1
              << matrix.value(k);
        }
      }

      // Go over non-core coefficients and compare them to those in the core.
      for (int k = matrix.col_start(i),
          end = matrix.col_start(i + 1); k != end; ++k) {
        int con_index = matrix.row_index(k);
        int scenario_index = con_info[con_index].scenario_index;
        if (scenario_index == 0)
          continue;
        int core_con_index = con_info[con_index].core_index;
        double coef = matrix.value(k);
        double core_coef = core_coefs[core_con_index];
        if (coef != core_coef) {
          int stage = stage_suffix.int_value(i) - 1;
          if (stage < 1) {
            // TODO: error: coefficient is inconsistent between scenarios
          }
          scenarios[scenario_index].Add(
              core_con_index, core_var_index, coef);
          if (core_coef == 0) {
            writer.Write("    C{:<7}  R{:<7}  0\n")
                << core_var_index + 1 << core_con_index + 1;
          }
        }
      }
    }

    writer.Write("RHS\n");
    for (int i = 0; i < num_cons; ++i) {
      if (con_info[i].scenario_index != 0)
        continue;
      double lb = p.con_lb(i), ub = p.con_ub(i);
      // TODO: ranges
      writer.Write("    RHS1      R{:<7}  {}\n")
          << con_info[i].core_index + 1 << (ub >= Infinity ? lb : ub);
    }

    writer.Write("BOUNDS\n");
    for (int i = 0; i < num_vars; ++i) {
      if (var_info[i].scenario_index != 0)
        continue;
      double lb = p.var_lb(i), ub = p.var_ub(i);
      if (lb != 0) {
        writer.Write(" LO BOUND1      C{:<7}  {}\n")
            << var_info[i].core_index + 1 << lb;
      }
      if (ub < Infinity) {
        writer.Write(" UP BOUND1      C{:<7}  {}\n")
            << var_info[i].core_index + 1 << ub;
      }
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
      writer.Write(" SC SCEN1     'ROOT'    {:<12}   T1\n") << probabilities[0];
      for (size_t i = 1, n = scenarios.size(); i < n; ++i) {
        writer.Write(" SC SCEN{:<4}  SCEN1     {:<12}   T2\n")
            << i + 1 << probabilities[i];
        for (auto t = scenarios[i].begin(),
            end = scenarios[i].end(); t != end; ++t) {
          writer.Write("    C{:<7}  R{:<7}  {}\n")
              << t->var_index + 1 << t->con_index + 1 << t->coef;
        }
        // TODO: write random rhs and bounds
      }
    }
    writer.Write("ENDATA\n");
  }
}
}
