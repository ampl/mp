/*
 AMPL to IBM/ILOG Concert interface.

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

#ifndef MP_SOLVERS_ILOGCP_CONCERT_H_
#define MP_SOLVERS_ILOGCP_CONCERT_H_

#ifdef __APPLE__
#include <limits.h>
#include <string.h>
#endif

#include <ilconcert/ilomodel.h>

#include <algorithm>  // std::find_if
#include <map>
#include <memory>

#ifdef MP_USE_HASH
# include <unordered_map>
#endif

#include "mp/expr-visitor.h"

namespace mp {

template <typename Alloc>
class BasicProblem;

class LinearExpr;

class EqualNumberOfArgs {
 public:
  bool operator()(IteratedExpr lhs, IteratedExpr rhs) const;
};

template <typename NumberOf>
class MatchNumberOfArgs {
 private:
  IteratedExpr impl_;

 public:
  explicit MatchNumberOfArgs(IteratedExpr e) : impl_(e) {}

  // Returns true if the stored expression has the same arguments as the nof's
  // expression.
  bool operator()(const NumberOf &nof) const {
    return EqualNumberOfArgs()(impl_, nof.expr);
  }
};

#ifdef MP_USE_HASH
class HashNumberOfArgs {
 public:
  std::size_t operator()(IteratedExpr e) const;
};
#endif

// A map from numberof expressions with the same argument lists to
// values and corresponding variables.
template <typename Var, typename CreateVar>
class NumberOfMap {
 public:
  typedef std::map<double, Var> ValueMap;

  struct NumberOf {
    IteratedExpr expr;  // numberof expression
    ValueMap values;

    explicit NumberOf(IteratedExpr e) : expr(e) {}
  };

 private:
  CreateVar create_var_;

#ifdef MP_USE_HASH
  // Map from a numberof expression to an index in numberofs_.
  typedef std::unordered_map<IteratedExpr,
    std::size_t, HashNumberOfArgs, EqualNumberOfArgs> Map;
  Map map_;
#endif

  std::vector<NumberOf> numberofs_;

 public:
  explicit NumberOfMap(CreateVar cv) : create_var_(cv) {}

  typedef typename std::vector<NumberOf>::const_iterator iterator;

  iterator begin() const {
    return numberofs_.begin();
  }

  iterator end() const {
    return numberofs_.end();
  }

  // Adds a numberof expression with a constant value.
  Var Add(double value, IteratedExpr e);
};

template <typename Var, typename CreateVar>
Var NumberOfMap<Var, CreateVar>::Add(double value, IteratedExpr e) {
  assert(Cast<NumericConstant>(e.arg(0)).value() == value);
#ifdef MP_USE_HASH
  std::pair<typename Map::iterator, bool> result =
      map_.insert(typename Map::value_type(e, numberofs_.size()));
  if (result.second)
    numberofs_.push_back(NumberOf(e));
  ValueMap &values = numberofs_[result.first->second].values;
# else
  typename std::vector<NumberOf>::reverse_iterator np =
      std::find_if(numberofs_.rbegin(), numberofs_.rend(),
                   MatchNumberOfArgs<NumberOf>(e));
  if (np == numberofs_.rend()) {
    numberofs_.push_back(NumberOf(e));
    np = numberofs_.rbegin();
  }
  ValueMap &values = np->values;
#endif  // MP_USE_HASH
  typename ValueMap::iterator i = values.lower_bound(value);
  if (i != values.end() && !values.key_comp()(value, i->first))
    return i->second;
  Var var(create_var_());
  values.insert(i, typename ValueMap::value_type(value, var));
  return var;
}

// Converter of optimization problems from NL to Concert format.
class MPToConcertConverter : public ExprVisitor<MPToConcertConverter, IloExpr> {
 private:
  IloEnv env_;
  IloModel model_;
  IloNumVarArray vars_;
  IloRangeArray cons_;
  std::vector<IloExpr> common_exprs_;
  unsigned flags_;

  class CreateVar {
   private:
    IloEnv env_;

   public:
    CreateVar(IloEnv env) : env_(env) {}

    IloIntVar operator()() const {
      return IloIntVar(env_, IloIntMin, IloIntMax);
    }
  };

  typedef NumberOfMap<IloIntVar, CreateVar> IlogNumberOfMap;
  IlogNumberOfMap numberofs_;

  IloNumExprArray ConvertArgs(VarArgExpr e);
  IloIntVar ConvertArg(CallExpr call, int index,
                       IloInt lb = IloIntMin, IloInt ub = IloIntMax);

  bool ConvertGlobalConstraint(CallExpr expr, IloConstraint &con);

  // Converts a sum of nonlinear and linear expressions into Concert format.
  IloExpr ConvertExpr(const LinearExpr &linear, NumericExpr nonlinear);

  // Converts a pairwise expression (alldiff or !alldiff).
  template <typename Constraint, bool negate>
  IloConstraint Convert(PairwiseExpr e) {
    int n = e.num_args();
    std::vector<IloExpr> args(n);
    int index = 0;
    for (PairwiseExpr::iterator i = e.begin(), end = e.end(); i != end; ++i)
      args[index++] = Visit(*i);
    Constraint alldiff(env_);
    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j)
        alldiff.add(negate ? args[i] == args[j] : args[i] != args[j]);
    }
    return alldiff;
  }

  // Converts logical expressions from NL to Concert format.
  class LogicalExprConverter :
      public ExprConverter<LogicalExprConverter, IloConstraint> {
   private:
    MPToConcertConverter &converter_;  // Main converter.

   public:
    explicit LogicalExprConverter(MPToConcertConverter &c) : converter_(c) {}

    using ExprConverter<LogicalExprConverter, IloConstraint>::Visit;

    IloExpr Visit(NumericExpr e) {
      return converter_.Visit(e);
    }

    IloConstraint VisitLogicalConstant(LogicalConstant c) {
      return IloNumVar(converter_.env_, 1, 1) == c.value();
    }

    IloConstraint VisitLT(RelationalExpr e) {
      return Visit(e.lhs()) < Visit(e.rhs());
    }

    IloConstraint VisitLE(RelationalExpr e) {
      return Visit(e.lhs()) <= Visit(e.rhs());
    }

    IloConstraint VisitEQ(RelationalExpr e) {
      return Visit(e.lhs()) == Visit(e.rhs());
    }

    IloConstraint VisitGE(RelationalExpr e) {
      return Visit(e.lhs()) >= Visit(e.rhs());
    }

    IloConstraint VisitGT(RelationalExpr e) {
      return Visit(e.lhs()) > Visit(e.rhs());
    }

    IloConstraint VisitNE(RelationalExpr e) {
      return Visit(e.lhs()) != Visit(e.rhs());
    }

    IloConstraint VisitOr(BinaryLogicalExpr e) {
      return IloIfThen(converter_.env_, !Visit(e.lhs()), Visit(e.rhs()));
    }

    IloConstraint VisitExists(IteratedLogicalExpr e);

    IloConstraint VisitAnd(BinaryLogicalExpr e) {
      return Visit(e.lhs()) && Visit(e.rhs());
    }

    IloConstraint VisitForAll(IteratedLogicalExpr e);

    IloConstraint VisitNot(NotExpr e) {
      return !Visit(e.arg());
    }

    IloConstraint VisitIff(BinaryLogicalExpr e) {
      return Visit(e.lhs()) == Visit(e.rhs());
    }

    IloConstraint VisitImplication(ImplicationExpr e);

    IloConstraint VisitAllDiff(PairwiseExpr e) {
      return converter_.Convert<IloAnd, false>(e);
    }

    IloConstraint VisitNotAllDiff(PairwiseExpr e) {
      return converter_.Convert<IloOr, true>(e);
    }
  };

 public:
  // Flags.
  enum {
    USENUMBEROF = 1,
    DEBUG       = 2
  };
  MPToConcertConverter(IloEnv env, unsigned flags);

  IloModel model() const { return model_; }
  IloNumVarArray vars() const { return vars_; }
  IloRangeArray cons() const { return cons_; }

  IloExpr Visit(NumericExpr e) {
    if ((flags_ & DEBUG) != 0)
      fmt::print("{}\n", str(e.kind()));
    return ExprVisitor<MPToConcertConverter, IloExpr>::Visit(e);
  }

  IloConstraint Visit(LogicalExpr e) {
    if ((flags_ & DEBUG) != 0)
      fmt::print("{}\n", str(e.kind()));
    return LogicalExprConverter(*this).Visit(e);
  }

  IloExpr VisitAdd(BinaryExpr e) {
    return Visit(e.lhs()) + Visit(e.rhs());
  }

  IloExpr VisitSub(BinaryExpr e) {
    return Visit(e.lhs()) - Visit(e.rhs());
  }

  IloExpr VisitMul(BinaryExpr e) {
    return Visit(e.lhs()) * Visit(e.rhs());
  }

  IloExpr VisitDiv(BinaryExpr e) {
    return Visit(e.lhs()) / Visit(e.rhs());
  }

  IloExpr VisitMod(BinaryExpr e) {
    IloNumExpr lhs(Visit(e.lhs())), rhs(Visit(e.rhs()));
    return lhs - IloTrunc(lhs / rhs) * rhs;
  }

  IloExpr VisitPow(BinaryExpr e) {
    return IloPower(Visit(e.lhs()), Visit(e.rhs()));
  }

  IloExpr VisitLess(BinaryExpr e) {
    return IloMax(Visit(e.lhs()) - Visit(e.rhs()), 0.0);
  }

  IloExpr VisitMin(VarArgExpr e) {
    return IloMin(ConvertArgs(e));
  }

  IloExpr VisitMax(VarArgExpr e) {
    return IloMax(ConvertArgs(e));
  }

  IloExpr VisitMinus(UnaryExpr e) {
    return -Visit(e.arg());
  }

  IloExpr VisitAbs(UnaryExpr e) {
    return IloAbs(Visit(e.arg()));
  }

  IloExpr VisitFloor(UnaryExpr e) {
    return IloFloor(Visit(e.arg()));
  }

  IloExpr VisitCeil(UnaryExpr e) {
    return IloCeil(Visit(e.arg()));
  }

  IloExpr VisitIf(IfExpr e);

  IloExpr VisitTanh(UnaryExpr e) {
    IloNumExpr exp(IloExponent(2 * Visit(e.arg())));
    return (exp - 1) / (exp + 1);
  }

  IloExpr VisitTan(UnaryExpr e) {
    return IloTan(Visit(e.arg()));
  }

  IloExpr VisitSqrt(UnaryExpr e) {
    return IloPower(Visit(e.arg()), 0.5);
  }

  IloExpr VisitSinh(UnaryExpr e) {
    IloNumExpr arg(Visit(e.arg()));
    return (IloExponent(arg) - IloExponent(-arg)) * 0.5;
  }

  IloExpr VisitSin(UnaryExpr e) {
    return IloSin(Visit(e.arg()));
  }

  IloExpr VisitLog10(UnaryExpr e) {
    return IloLog10(Visit(e.arg()));
  }

  IloExpr VisitLog(UnaryExpr e) {
    return IloLog(Visit(e.arg()));
  }

  IloExpr VisitExp(UnaryExpr e) {
    return IloExponent(Visit(e.arg()));
  }

  IloExpr VisitCosh(UnaryExpr e) {
    IloNumExpr arg(Visit(e.arg()));
    return (IloExponent(arg) + IloExponent(-arg)) * 0.5;
  }

  IloExpr VisitCos(UnaryExpr e) {
    return IloCos(Visit(e.arg()));
  }

  IloExpr VisitAtanh(UnaryExpr e) {
    IloNumExpr arg(Visit(e.arg()));
    return IloLog((1 + arg) / (1 - arg)) * 0.5;
  }

  IloExpr VisitAtan2(BinaryExpr e);

  IloExpr VisitAtan(UnaryExpr e) {
    return IloArcTan(Visit(e.arg()));
  }

  IloExpr VisitAsinh(UnaryExpr e) {
    IloNumExpr arg(Visit(e.arg()));
    return IloLog(arg + IloPower(IloSquare(arg) + 1, 0.5));
  }

  IloExpr VisitAsin(UnaryExpr e) {
    return IloArcSin(Visit(e.arg()));
  }

  IloExpr VisitAcosh(UnaryExpr e) {
    IloNumExpr arg(Visit(e.arg()));
    return IloLog(arg + IloPower(arg + 1, 0.5) * IloPower(arg - 1, 0.5));
  }

  IloExpr VisitAcos(UnaryExpr e) {
    return IloArcCos(Visit(e.arg()));
  }

  IloExpr VisitSum(SumExpr e);

  IloExpr VisitTruncDiv(BinaryExpr e) {
    return IloTrunc(Visit(e.lhs()) / Visit(e.rhs()));
  }

  IloExpr VisitRound(BinaryExpr e);

  IloExpr VisitTrunc(BinaryExpr e);

  IloExpr VisitCount(CountExpr e);

  IloExpr VisitNumberOf(NumberOfExpr e);

  IloExpr VisitPowConstExp(BinaryExpr e) {
    return IloPower(Visit(e.lhs()), Cast<NumericConstant>(e.rhs()).value());
  }

  IloExpr VisitPow2(UnaryExpr e) {
    return IloSquare(Visit(e.arg()));
  }

  IloExpr VisitPowConstBase(BinaryExpr e) {
    return IloPower(Cast<NumericConstant>(e.lhs()).value(), Visit(e.rhs()));
  }

  IloExpr VisitPLTerm(PLTerm e);

  IloExpr VisitCall(CallExpr e);

  IloExpr VisitNumericConstant(NumericConstant n) {
    return IloExpr(env_, n.value());
  }

  IloExpr VisitVariable(Reference r) {
    return vars_[r.index()];
  }

  IloExpr VisitCommonExpr(Reference r) {
    return common_exprs_[r.index()];
  }

  // Combines 'numberof' operators into IloDistribute constraints
  // which are much more useful to the solution procedure.
  void FinishBuildingNumberOf();

  typedef BasicProblem< std::allocator<char> > Problem;

  void Convert(const Problem &p);

  /// [[ The interface ]]
  void PushVariables(const Problem& p);
  void PushCommonSubExpr(const Problem& p);
  void PushObjectives(const Problem& p);
  void PushAlgebraicConstraints(const Problem& p);
  void PushLogicalConstraints(const Problem& p);
  void FinishConversion(const Problem& p);

};
}

#endif  // MP_SOLVERS_ILOGCP_CONCERT_H_
