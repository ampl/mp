#include "baronmpmodelapi.h"

#include <limits>
namespace mp {

  std::function<std::string(int)> VExpr::varName = nullptr;

 

  void BaronmpModelAPI::InitProblemModificationPhase(
    const FlatModelInfo*) {

    // Have to be written after eventual additional options are injected
    initBaronFile();
    writeBaronOptions();
    VExpr::varName = std::bind(&BaronmpModelAPI::varName, this, std::placeholders::_1);
  }


  std::string sanitizeName(fmt::CStringRef name , char prefix) {
    std::string n(name.c_str());
    // Check if the first character is a letter; if not, prepend the prefix
    if (!std::isalpha(n[0])) {
      n = prefix + n;
    }

    // Replace non-alphanumeric characters or spaces with underscores
    for (char& ch : n) {
      if (!std::isalnum(ch)) {
        ch = '_';
      }
    }
    return n;
  }


  void BaronmpModelAPI::AddBounds(const double* bounds, 
    const std::vector<VTYPE>& vtypes,bool lower) {
    bool headerWritten = false;

    double bnd = lower ? -std::numeric_limits<double>::infinity() :
      std::numeric_limits<double>::infinity();
    
    for (int i = 0; i < vtypes.size(); i++)
    {
      if (vtypes[i] == VTYPE::BIN)
        continue;
      if (lower &&   bounds[i]==0 && vtypes[i] == VTYPE::POS)
        continue;
      if (bounds[i] != bnd)
      {
        if (!headerWritten) {
          writeVars(vars_buffer, lower ? "LOWER BOUNDS\{\n" : "UPPER_BOUNDS\{\n");
          headerWritten = true;
        }
        writeVars(vars_buffer, fmt::format("{}: {:.{}g};\n", lp()->varNames[i], bounds[i], std::numeric_limits<double>::max_digits10));
      }
    }
    if (headerWritten)
      writeVars(vars_buffer, "}\n\n");
  }
  void BaronmpModelAPI::AddVariables(const std::vector<int>& indices, 
    fmt::StringRef prefix, fmt::StringRef header, const char* const* names) {
    
    if (indices.size() > 0) {
      std::string name;
      writeVars(vars_buffer, fmt::format("{} ", header.to_string()));
      for (int i = 0; i < indices.size(); i++)
      {
        if (i != 0)
          writeVars(vars_buffer, ", ");
        if (names)

          name = sanitizeName(names[indices[i]], 'v');
        else
          name = fmt::format("{}{}", prefix, indices[i]);
        lp()->varMap[name] = indices[i];
        lp()->varNames[indices[i]] = name;
        lp()->baronToAMPLIndices.push_back(indices[i]);

        writeVars(vars_buffer, name);
      }
      writeVars(vars_buffer, ";\n");
    }
  }

  void BaronmpModelAPI::AddVariables(const VarArrayDef& v) {
    lp()->varNames.resize(v.size());
    std::vector<int> indicesFree, indicesPos, indicesBin, indicesInt;
    std::vector<VTYPE> vtype(v.size());
    for (auto i = 0; i < v.size(); i++)
    {
      if (v.ptype()[i] == var::Type::INTEGER)
      {
        if ((v.pub()[i] == 1) && (v.plb()[i] == 0))
        {
          vtype[i] = VTYPE::BIN; // BIN
          indicesBin.push_back(i);
        }
        else
        {
          vtype[i] = VTYPE::INT; // INT
          indicesInt.push_back(i);
        }
      }
      else
      {
        if (v.plb()[i] == 0)
        {
          vtype[i] = VTYPE::POS; // POS
          indicesPos.push_back(i);
        }
        else
        {
          vtype[i] = VTYPE::FREE; // FREE
          indicesFree.push_back(i);
        }
      }
    
    }
    int cv = 0;
    std::string name;
    AddVariables(indicesBin, "b", "BINARY_VARIABLES", v.pnames());
    AddVariables(indicesInt, "i", "INTEGER_VARIABLES", v.pnames());
    AddVariables(indicesPos, "p", "POSITIVE_VARIABLES", v.pnames());
    AddVariables(indicesFree, "x", "VARIABLES", v.pnames());
    AddBounds(v.plb(), vtype, true);
    AddBounds(v.pub(), vtype, false);
    nVarsBinary(indicesBin.size());
    nVarsContinuous(indicesPos.size() + indicesFree.size());
    nVarsInteger(indicesInt.size());
  }

  void addCoefficient(fmt::MemoryWriter &w, double coeff, bool first)
  {
    if (coeff > 0) {
        if (!first) w << " + ";
        // Omit '1*' for positive coefficient 1
        if (coeff != 1) 
            w << coeff << "*";
    } 
    else { // Coefficient is negative
        // If the coefficient is -1, just print ' - ', otherwise print the coefficient
        if (coeff == -1)  w << " - ";
         else  w << coeff << "*";
    }
  }

  void appendLinearTerm(fmt::MemoryWriter &w, double coeff, const char* vname, bool first = false)
  {
    addCoefficient(w, coeff, first);
    w << vname;
  }
  void appendQuadTerm(fmt::MemoryWriter &w, double coeff, const char* v1, const char* v2=nullptr, bool first = false)
  {
    addCoefficient(w, coeff, first);
    if(v2)
      w << v1 << "*" << v2;
    else
      w << v1 << "^2";
  }


  void BaronmpModelAPI::AddObjectiveHeader(fmt::MemoryWriter& w, int sense) {
    w << "OBJ: ";
    if (sense == mp::obj::Type::MAX)
      w << "maximize ";
    else
      w << "minimize ";
    lp()->obj_sense = sense;

  }
  void BaronmpModelAPI::SetLinearObjective(int iobj, const LinearObjective& lo) {
    if (iobj < 1) {
      fmt::MemoryWriter w;
      AddObjectiveHeader(w, lo.obj_sense());
      
      auto nvars = lo.GetLinTerms().size();
      auto coeffs = lo.GetLinTerms().pcoefs();
      auto vars = lo.GetLinTerms().pvars();

      for (size_t i = 0; i < nvars; i++)
        appendLinearTerm(w, coeffs[i], varName(vars[i]).c_str(), i == 0);
      w<<(";\n");
      obj = w.str();
     
    }
    else {
      MP_RAISE("Baron supports only one objective.");
    }
  }

  void BaronmpModelAPI::SetQuadraticObjective(int iobj, const QuadraticObjective& qo) {
    if (iobj < 1) {

      auto lt = qo.GetLinTerms();
      auto qt = qo.GetQPTerms();
      fmt::MemoryWriter w;
      AddObjectiveHeader(w, qo.obj_sense());

      size_t i;
      for (i = 0; i < lt.size(); i++)
        appendLinearTerm(w, lt.pcoefs()[i], varName(lt.var(i)).c_str(), i == 0);
      bool addedLinearTerms = i > 0;
      for (i = 0; i < qt.size(); i++)
        appendQuadTerm(w, qt.coef(i), varName(qt.var1(i)).c_str(),
          (qt.var1(i) == qt.var2(i)) ? nullptr : varName(qt.var2(i)).c_str(), (!addedLinearTerms) && (i == 0));
      w << ";\n";
      obj = w.str();
    }
    else {
      MP_RAISE("Baron supports only one objective.");
    }
  }
  void BaronmpModelAPI::SetNLObjective(int i , const NLObjective& nlo) {
    if(i > 0) MP_RAISE("Baron supports only one objective.");

    const auto& exp = GetExpression(nlo);
    fmt::MemoryWriter w;
    AddObjectiveHeader(w, nlo.obj_sense());

    auto lin = nlo.GetLinTerms();
    
    for (int i = 0; i < lin.size(); ++i)
      appendLinearTerm(w, lin.coef(i), varName(lin.var(i)).c_str(), i == 0);
    if (lin.size() > 0) w << " + ";
    exp.append(w, false);
    w << ";\n";
    obj = w.str();
  }

  std::string BaronmpModelAPI::varName(int index) {
    return lp()->varNames[index];
  }

  std::string BaronmpModelAPI::createConName(const std::string& name) {
    std::string n;
    if (name.empty())
      n = fmt::format("c{}", ++n_unamed_constraint);
    else
      n = sanitizeName(name, 'c');
    lp()->conNames.push_back(n);
    return n;
  }

  void appendNameAndLhs(fmt::MemoryWriter &w, const std::string &name, double lhs, double rhs) {
    w << name << ": ";
    bool hasLhs = lhs != -std::numeric_limits<double>::infinity();
    bool hasRhs = rhs != std::numeric_limits<double>::infinity();
    if ((hasLhs && hasRhs) && (lhs != rhs))
      w << lhs << " <= ";

  }
  void appendRhs(fmt::MemoryWriter &w, double lhs, double rhs) {
     bool hasLhs = lhs != -std::numeric_limits<double>::infinity();
    bool hasRhs = rhs != std::numeric_limits<double>::infinity();
    if (hasRhs)
      w <<  ((lhs != rhs) ? "<=" : "==") << rhs;
    else if (hasLhs) // range
      w << " >= " << lhs;
    w << ";\n";

  }

  void BaronmpModelAPI::addLinear(const std::string& name,
    double lhs, double rhs,
    size_t nvars, const int* vars, const double* coeffs){
    if (nvars == 0)
      return;
    fmt::MemoryWriter w;
    appendNameAndLhs(w, createConName(name), lhs, rhs);
    for (size_t i = 0; i < nvars; i++)
      appendLinearTerm(w, coeffs[i], varName(vars[i]).c_str(), i == 0);
    appendRhs(w, lhs, rhs);
    cons.push_back(w.str());
}

  void BaronmpModelAPI::addQuadratic(const std::string& name,
    double lhs, double rhs, const mp::LinTerms& lt,
    const mp::QuadTerms& qt) {

    bool addedLinearTerm = false;
    fmt::MemoryWriter w;
    appendNameAndLhs(w, createConName(name), lhs, rhs);
    size_t i;
    for (i = 0; i < lt.size(); i++) 
      appendLinearTerm(w, lt.pcoefs()[i], varName(lt.var(i)).c_str(), i == 0);
    bool addedLinearTerms = i > 0;
    for (i = 0; i < qt.size(); i++)
       appendQuadTerm(w, qt.coef(i), varName(qt.var1(i)).c_str(), 
          (qt.var1(i) == qt.var2(i)) ? nullptr : varName(qt.var2(i)).c_str(), (!addedLinearTerms) && (i==0));
    appendRhs(w, lhs, rhs);
    cons.push_back(w.str());
  }

  template <typename Args, typename Params, typename NumOrLogic, typename Id>
  void BaronmpModelAPI::addFunctionalConstraint(const std::string& func,
    const CustomFunctionalConstraint<Args, Params, NumOrLogic, Id>& c) {
    fmt::MemoryWriter w;
    
    w<< createConName(c.GetName()) << ": ";

    w << fmt::format("{} - {}({}", varName(c.GetResultVar()), func, varName(c.GetArguments()[0]));
    if (c.GetArguments().size() > 1) {
      for (auto i = 1; i < c.GetArguments().size(); ++i)
        w << fmt::format(", {}", varName(c.GetArguments()[i]));
    }
    w << ") = 0;\n";
    cons.push_back(w.str());

  }


void BaronmpModelAPI::AddConstraint(const LinConRange& lc) {
  addLinear(lc.name(), lc.lb(), lc.ub(),
    lc.size(), lc.pvars(), lc.pcoefs()) ;
}
void BaronmpModelAPI::AddConstraint(const LinConLE& lc) {
  addLinear(lc.name(), 
    -std::numeric_limits<double>::max(),lc.ub(),
    lc.size(), lc.pvars(), lc.pcoefs());
}
void BaronmpModelAPI::AddConstraint(const LinConEQ& lc) {
  addLinear(lc.name(),
    lc.lb(), lc.ub(),
    lc.size(), lc.pvars(), lc.pcoefs());
}
void BaronmpModelAPI::AddConstraint(const LinConGE& lc) {
  addLinear(lc.name(),
    lc.lb(), std::numeric_limits<double>::max(),
    lc.size(), lc.pvars(), lc.pcoefs());
}

void BaronmpModelAPI::AddConstraint(const QuadConRange& qc) {
  addQuadratic(qc.name(), qc.lb(), qc.ub(), qc.GetLinTerms(), qc.GetQPTerms());
}

void BaronmpModelAPI::AddConstraint( const QuadConLE& qc ) {
  addQuadratic(qc.name(), qc.lb(), qc.ub(), qc.GetLinTerms(), qc.GetQPTerms());
}

void BaronmpModelAPI::AddConstraint( const QuadConEQ& qc ) {
  addQuadratic(qc.name(), qc.lb(), qc.ub(), qc.GetLinTerms(), qc.GetQPTerms());
}

void BaronmpModelAPI::AddConstraint( const QuadConGE& qc ) {
  addQuadratic(qc.name(), qc.lb(), qc.ub(), qc.GetLinTerms(), qc.GetQPTerms());
}


void BaronmpModelAPI::AddConstraint(const ExpConstraint& cc) {
  addFunctionalConstraint("exp", cc);
}

void BaronmpModelAPI::AddConstraint(const LogConstraint& cc) {
  addFunctionalConstraint("log", cc);
}

void BaronmpModelAPI::AddConstraint(const ExpAConstraint& cc) {
  // a^v
  fmt::MemoryWriter w;
  w<< createConName(cc.GetName()) << ": ";
  w << fmt::format("{} = {} ^ {};\n",  varName(cc.GetResultVar()), cc.GetParameters()[0],varName( cc.GetArguments()[0]));
  cons.push_back(w.str());
}

void BaronmpModelAPI::AddConstraint(const LogAConstraint& cc) {
  fmt::MemoryWriter w;
  w<< createConName(cc.GetName()) << ": ";
  w << fmt::format("{} = log({})*;\n", varName(cc.GetResultVar()), 
    varName( cc.GetArguments()[0]), 1/std::log(cc.GetParameters()[0]));
  cons.push_back(w.str());
}

void BaronmpModelAPI::AddConstraint(const PowConstExpConstraint& cc) {
  // v ^ a
   fmt::MemoryWriter w;
  w<< createConName(cc.GetName()) << ": ";
  w << fmt::format("{} = {} ^ {};\n",  varName(cc.GetResultVar()),varName( cc.GetArguments()[0]), cc.GetParameters()[0]);
  cons.push_back(w.str());
}


template <int SENSE> void BaronmpModelAPI::addTopLevel(const NLBaseAssign<SENSE>& c) {
  const auto var = GetVariable(c);
  fmt::MemoryWriter w;
  w << createConName(c.GetName()) << ": ";
  const std::string sense[] = { "<=", "==", ">="};
  // Baron does not accept x >= f(x)
  //w << fmt::format("{} {} ", varName(var), sense[SENSE+1]);
  w << fmt::format("{} - ", varName(var));
  const auto& frm = GetExpression(c);
  w << frm.ToString(false);
  w << fmt::format("{} 0", sense[SENSE + 1]);
  

  w << ";\n";
  cons.push_back(w.str());
}

void BaronmpModelAPI::AddConstraint(const NLConstraint& nlc) {
  const auto& exp = GetExpression(nlc);
  double lhs=GetLower(nlc), rhs=GetUpper(nlc);
  fmt::MemoryWriter w;
  appendNameAndLhs(w, createConName(nlc.GetName()), lhs, rhs);
  for (int i = 0; i < GetLinSize(nlc); ++i)
      appendLinearTerm(w, GetLinCoef(nlc, i), varName(GetLinVar(nlc, i)).c_str(), i == 0);
  if (GetLinSize(nlc) > 0) w <<" + ";
  exp.append(w, false);
  appendRhs(w, lhs, rhs);
  cons.push_back(w.str());
}

void BaronmpModelAPI::AddConstraint(const NLAssignEQ& nla) {
  addTopLevel(nla);
}
void BaronmpModelAPI::AddConstraint(const NLAssignLE& nla) {
  addTopLevel(nla);
}
void BaronmpModelAPI::AddConstraint(const NLAssignGE& nla) {    
  addTopLevel(nla);
}

VExpr BaronmpModelAPI::AddExpression(const NLAffineExpression& nla) {
  VExpr tree;
  tree.opcode = Opcode::ADD;
  AppendLinAndConstTerms(tree, nla);
  return tree;
}

VExpr  BaronmpModelAPI::AddExpression(const NLQuadExpression& nlq) {
  VExpr tree;
  tree.opcode = Opcode::ADD;
  AppendLinAndConstTerms(tree, nlq);
  AppendQuadTerms(tree, nlq);
  return tree;
}

template <class MPExpr>
void BaronmpModelAPI::AppendLinAndConstTerms(
  Expr &ff, const MPExpr& nla) {
  if (double ct = GetConstTerm(nla))
    ff.AddArgument(VExpr::makeConstant(ct));

  for (int i = 0; i < GetLinSize(nla); ++i) {
    if (double coef = GetLinCoef(nla, i)) {
      if (1.0 != coef) {
        auto f1 = VExpr(Opcode::MULTIPLY);
        f1.AddArgument(VExpr::makeConstant(coef));
        f1.AddArgument(GetLinTerm(nla, i));
        ff.AddArgument(f1);
      }
      else {
        ff.AddArgument(GetLinTerm(nla, i));
      }
    }
  }
}


template <class MPExpr>
void BaronmpModelAPI::AppendQuadTerms(
  VExpr& ff, const MPExpr& nlq) {
  for (int i = 0; i < GetQuadSize(nlq); ++i) {
    if (double coef = GetQuadCoef(nlq, i)) {
      auto f1 = VExpr(Opcode::MULTIPLY);
      f1.AddArgument(GetQuadTerm1(nlq, i));
      f1.AddArgument(GetQuadTerm2(nlq, i));
      if (1.0 != coef) {
        f1.AddArgument(VExpr::makeConstant(coef));
      }
      ff.AddArgument(f1);
    }
  }
}

VExpr BaronmpModelAPI::AddExpression(const DivExpression& e) {
  auto res = VExpr(Opcode::DIVIDE);
  res.AddArgument(GetArgExpression(e, 0));
  res.AddArgument(GetArgExpression(e, 1));
  return res;

}
VExpr BaronmpModelAPI::AddExpression(const ExpAExpression& e) {
  auto res = VExpr(Opcode::POW);
  // Constant^Variable
  res.AddArgument(VExpr::makeConstant(GetParameter(e, 0)));
  res.AddArgument(GetArgExpression(e, 0));
  return res;
}

VExpr BaronmpModelAPI::AddExpression(const ExpExpression& e) {
  // e^v
  return VExpr(Opcode::EXP, GetArgExpression(e, 0));
}

VExpr BaronmpModelAPI::AddExpression(const LogExpression& e) {
  return VExpr(Opcode::LOG, GetArgExpression(e, 0));
}
VExpr BaronmpModelAPI::AddExpression(const LogAExpression& e) {
  auto logvar = VExpr(Opcode::LOG, GetArgExpression(e, 0));
  auto logbase = VExpr::makeConstant(1/std::log(GetParameter(e, 0)));
  auto expr = VExpr(Opcode::MULTIPLY);
  expr.AddArgument(logvar);
  expr.AddArgument(logbase);
  return expr;
}

VExpr BaronmpModelAPI::AddExpression(const PowConstExpExpression& e) {
  auto expr=VExpr(Opcode::POW, GetArgExpression(e, 0));
  for (int i = 0; i < GetNumParameters(e); ++i)
    expr.AddArgument(VExpr::makeConstant(GetParameter(e, i)));
  return expr;
}



void BaronmpModelAPI::FinishProblemModificationPhase() {
  // Write out vars
  writeBaron(vars_buffer.str());

  // Write out constraints
  if (lp()->conNames.size() != 0)
  {
    writeBaron("EQUATIONS ");
    writeBaron(lp()->conNames[0]);
    if(lp()->conNames.size() > 1)
    for (auto i = 1; i < lp()->conNames.size(); ++i)
      writeBaron(" ,{}", lp()->conNames[i]);
    writeBaron(";\n\n");

    for (auto c : cons)
      writeBaron(c);
  }
  writeBaron("\n");

  // Write objective
  if (obj.empty()) // Baron does not like no objectives
    writeBaron("OBJ:\nminimize 0; ");
  else
    writeBaron(obj);
}




} // namespace mp
