#ifndef MP_CUOPT_SOLVERMODEL_H_
#define MP_CUOPT_SOLVERMODEL_H_

#include <ostream>
#include <vector>
#include <variant>
#include <string>
#include <map>
#include <unordered_map>
#include <mp/format.h>

// This file is normally a header file coming from the solver library
// distribution
namespace Solver {


  enum ATTRIBS {
    NVARS, 
    NVARS_INT,
    NVARS_CONT,
    NVARS_BIN,

    NCONS,
    NCONS_TYPE,
    NEXPR_TYPE,

    NOBJS,
    ISQOBJ
  };

  
  /// Enumeration to differentiate flat constraint types
  /// Added with the ModelApi's addConstraint overloads
  enum ConsType {
    CONS_LIN,
    CONS_QUAD,
    CONS_QUAD_CONE,
    CONS_QUAD_CONE_ROTATED,
    CONS_QUAD_CONE_EXP,
    CONS_INDIC,
    CONS_SOS,

    CONS_MAX,
    CONS_MIN,
    CONS_ABS,
    CONS_AND,
    CONS_OR,

    CONS_EXP,
    CONS_EXPA,
    CONS_LOG,
    CONS_LOGA,

    CONS_POW,
    CONS_POW_VAR,
    CONS_SIN,
    CONS_COS,
    CONS_TAN,

    CONS_PL,
    CONS_NL
  };


  // Define a Node (VExpr) for the expression tree
  class VExpr {
  public:
    // Capturing expression trees 
    enum class Opcode {
      CONSTANT,
      VARIABLE,
      ADD,
      SUBTRACT,
      MULTIPLY,
      DIVIDE,

      

      EXP,
      POW,
      LOGA,
      ABS,
      LOG,


      SIN,
      COS,
      TAN,
      SINH,
      COSH,
      TANH,
      ASIN,
      ACOS,
      ATAN,
      ASINH,
      ACOSH,
      ATANH
    };
    static inline const std::unordered_map<Opcode, std::string> opcodeStrings = {
        {Opcode::CONSTANT, "CONSTANT"},
        {Opcode::VARIABLE, "VARIABLE"},
        {Opcode::ADD, "+"},
        {Opcode::SUBTRACT, "-"},
        {Opcode::MULTIPLY, "*"},
        {Opcode::DIVIDE, "/"},
        {Opcode::POW, "^"},
        {Opcode::LOGA, "loga"},
        {Opcode::ABS, "abs"},
        {Opcode::EXP, "exp"},
        {Opcode::LOG, "log"},
        {Opcode::SIN, "sin"},
        {Opcode::COS, "cos"},
        {Opcode::TAN, "tan"},
        {Opcode::SINH, "sinh"},
        {Opcode::COSH, "cosh"},
        {Opcode::TANH, "tanh"},
        {Opcode::ASIN, "asin"},
        {Opcode::ACOS, "acos"},
        {Opcode::ATAN, "atan"},
        {Opcode::ASINH, "asinh"},
        {Opcode::ACOSH, "acosh"},
        {Opcode::ATANH, "atanh"}
    };
    void reset() {
      children.clear();
    }
    Opcode opcode;  // Stores the type of the operation or expression
    std::variant<double, int> value;  // Either a constant (double) or a variable index (int)
    std::vector<VExpr> children;  // Can store multiple children for n-ary operations
    VExpr* parent = nullptr;

    VExpr() : opcode(Opcode::CONSTANT), value(-1) {}

    VExpr(Opcode op) : opcode(op), value(-1) {
      assert((op != Opcode::CONSTANT) && (op != Opcode::VARIABLE));
    }
    VExpr(Opcode op, VExpr arg) : opcode(op), value(-1) {
      AddArgument(arg);
    }
  private:
    VExpr(Opcode op, std::variant<double, int> val) : opcode(op), value(val) {}
  public:
    static VExpr makeVariable(int index) {
      return VExpr(Opcode::VARIABLE, index);
    }
    static VExpr makeConstant(double value) {
      return VExpr(Opcode::CONSTANT, value);
    }
    // Add an argument (child node) to this expression
    void AddArgument(VExpr child) {
      child.parent = this;
      children.push_back(child);
    }

  };

  enum ObjType {
    OBJ_LIN, 
    OBJ_QUAD,
    OBJ_NL
  };

  enum VarType {
    CONTINUOUS,
    INTEGER,
    BINARY
  };

  /// <summary>
  /// MP has 4 types of top level non-linear constraints:
  /// NLConstraint    lhs <= expr <= rhs
  /// NLAssignEQ      v = expr
  /// NLAssignGE      v >= expr
  /// NLAssignLE      v <= expr
  /// </summary>
  enum NLConsType {
    NLCons,
    AssignEQ,
    AssignGE,
    AssignLE,
    NLLogical,
    NLReifEquiv,
    NLReifImpl,
    NLReifRimpl
  };
  class SolverModel {
    
  public:
    struct Variable {
      std::string name;
      VarType type;
      double lb, ub;
    };

    struct Constraint {
      std::string name;
      ConsType type;
      double lhs, rhs;
    };
    struct ConstraintFlat : public Constraint {

    };

    struct ConstraintNL  : public Constraint{
      VExpr e;
      NLConsType nl_type;
      int resultvar; // used only for Assign NL constraints
      
    };

    struct Objective {
      std::string name;
      ObjType type;
    };
    struct ObjectiveFlat : public Objective {
      
    };

    struct ObjectiveNL : public Objective {
      VExpr e;
    };

    /// Used to accumulate statistics on the model
    /// Also used when creating default names for variables
    int nVarsPerType[3];
    const char prefix[3] = { 'x', 'i', 'b' };
    std::string create_conname(std::string_view name = "")
    {
      if (name.empty()) {
        return fmt::format("c{}", cons_.size());
      }
      return std::string(name); // Convert string_view to string if provided.
    }
    std::string create_varname(VarType t, std::string_view name = "")
    {
      if (name.empty()) {
        return fmt::format("{}{}", prefix[t], nVarsPerType[t]);
      }
      return std::string(name); // Convert string_view to string if provided.
    }
    std::string create_objname(std::string_view name = "")
    {
      if (name.empty()) {
        return fmt::format("o{}", objs_.size());
      }
      return std::string(name); // Convert string_view to string if provided.
    }
    std::map<ConsType, int> nEntities_;
    std::vector<Variable> vars_;
    std::vector< Constraint> cons_;
    std::vector<Objective> objs_;
    int verbosity_ = 1;

  public:
    void SetVerbosity(int level) { verbosity_ = level; }
    int GetVerbosity() { return verbosity_;  }
    int AddVariable(int type, double lb, double ub, std::string_view name = "") {
      Variable v;
      // Convert mp variable type to solver.
      // Sometimes solvers required explicit declaration for binary variables, 
      // so we detect binaries by finding integer vars with bounds
      // 0<=x<=1, 0<=x<=0, 1<=x<=1
      if (type == 0)
        v.type = VarType::CONTINUOUS;
      else {
        if (((lb == 0) && (ub == 1)) || ((lb == 0) && (ub == 0)) || ((lb == 1) && (ub == 1)))
          v.type = VarType::BINARY;
        else
          v.type = VarType::INTEGER;
      }
      nVarsPerType[v.type]++; // Accumulate number of variable
      v.name = create_varname(v.type, name);
      v.lb = lb,
      v.ub = ub;
      vars_.push_back(v);
      return 0;
    }

    std::string AddConstraintFlat(ConsType t, std::string_view name = "") {
      ConstraintFlat c = { create_conname(name), t, 0,0};
      cons_.push_back(c);
      if (nEntities_.find(t) != nEntities_.end())
        nEntities_[t] += 1;
      else
        nEntities_[t] = 1;
      return c.name;
    }
    std::string AddConstraintNLAssign(int resvar, NLConsType nl_type, VExpr e, std::string_view name = "") {
      ConstraintNL c = { create_conname(name), ConsType::CONS_NL, 0, 0, e , nl_type, resvar};
      cons_.push_back(c);
      fmt::MemoryWriter w;
      const std::map<NLConsType, std::string> op = { {AssignEQ, "=="}, {AssignGE, ">="}, {AssignLE, "<="}, {NLReifEquiv, "<==>"},
        {NLReifImpl, "==>"}, {NLReifRimpl, "<=="}};
      w << c.name << ": ";

      if (nl_type == NLLogical) // NLLogical does not need a variable 
        w << "true = ";
      else {
        // Append variable name
        w << var_name(resvar).data();
        // For implications, we show the comparison explicitly
        if ((nl_type == NLReifImpl) || (nl_type == NLReifRimpl))
          w << "==1 ";
        // Append operand
        w << " " << op.at(nl_type) << " ";
      }
      // Append expression
      append(w, e, true);
      fmt::print(w.str());
      return c.name;
    }
    std::string AddConstraintNL(VExpr e, double lhs, double rhs, std::string_view name = "") {
      ConstraintNL c = { create_conname(name), ConsType::CONS_NL, rhs, lhs, e, NLCons };
      cons_.push_back(c);
      fmt::MemoryWriter w;
      w << c.name << ": " << lhs << " <= ";
      append(w, e, false);
      w << " <= " << rhs << "\n";
      fmt::print(w.str());
      return c.name;
    }
    std::string AddObjective(std::string_view name, ObjType t) {
      ObjectiveFlat o = { create_objname(name), t };
      objs_.push_back(o);
      return o.name;
    }
    std::string AddNLObjective(std::string_view name, VExpr v, bool maximize) {
      ObjectiveNL o = { create_objname(name), ObjType::OBJ_NL, v };
      objs_.push_back(o);
      fmt::MemoryWriter w;
      w << o.name << " " << (maximize ? "maximize " : "minimize ");
      append(w, v, false);
      w << "\n";
      fmt::print(w.str());
      return o.name;
    }

    std::string_view var_name(std::size_t index) {
      return vars_[index].name;
    }
    int countObj(ObjType ot) {
      int i = 0;
      for (const auto& o : objs_)
        if (o.type == ot) i++;
      return i;


    }
    int GetAttribute(ATTRIBS a, ConsType t = ConsType::CONS_ABS) {
      switch (a) {
      case NCONS:
        return cons_.size();
      case NCONS_TYPE:
        return nEntities_[t];
      case NVARS:
        return vars_.size();
      case NVARS_CONT:
        return getNumVars(VarType::CONTINUOUS);
      case NVARS_INT:
        return getNumVars(VarType::INTEGER);
      case NVARS_BIN:
        return getNumVars(VarType::BINARY);
      case NOBJS:
        if(t == ConsType::CONS_ABS) // no value specified, total number of objs
          return objs_.size();
        if (t == ConsType::CONS_LIN) return countObj(ObjType::OBJ_LIN);
        if (t == ConsType::CONS_QUAD) return countObj(ObjType::OBJ_QUAD);
        if (t == ConsType::CONS_NL) return countObj(ObjType::OBJ_NL);
        throw std::runtime_error("Unexpected type when retrieving objective statistics");
      }
      return -1;

    }
    
    void allocateVars(int nvars) {
      vars_.reserve(nvars);
    }

    int getNumVars(VarType t) {
      return nVarsPerType[t];
    }


    void append(fmt::MemoryWriter& w, VExpr exp,
      bool endl = true) const {
      VExpr::Opcode opcode = exp.opcode;
      assert(opcode >= VExpr::Opcode::CONSTANT && opcode <= VExpr::Opcode::ATANH);
      if (opcode >= VExpr::Opcode::EXP) // 1-nary functions, print its name
        w << VExpr::opcodeStrings.at(opcode);
      if (!exp.children.empty()) { // print operands
        w << "(";
        for (size_t i = 0; i < exp.children.size(); ++i) {
          append(w, exp.children[i]);
          if (i < exp.children.size() - 1) {
            switch (opcode) {
            case VExpr::Opcode::ADD: w << " + "; break;
            case VExpr::Opcode::SUBTRACT: w << " - "; break;
            case VExpr::Opcode::MULTIPLY: w << " * "; break;
            case VExpr::Opcode::DIVIDE: w << " / "; break;
            case VExpr::Opcode::POW: w << "^"; break;
            case VExpr::Opcode::LOGA: w << ","; break;
            default: break;
            }
          }
        }
        w << ")";
      }
      else {
        if (opcode == VExpr::Opcode::CONSTANT) {
          w << std::get<double>(exp.value);
        }
        else if (opcode == VExpr::Opcode::VARIABLE) {
          w << vars_[std::get<int>(exp.value)].name;
        }

      }
      if (endl && (exp.parent == nullptr))
        w << "\n";
    }

    // Print the expression in-order
    std::string ToString(VExpr exp,bool endl = true) const {
      fmt::MemoryWriter w;
      append(w, exp, endl);
      return w.str();
    }


   
  };

  SolverModel* CreateSolverModel();

}

#endif // MP_CUOPT_SOLVERMODEL_H_
