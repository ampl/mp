/**
 * Convert various constraints to cones.
 *
 * These conversions might be run before and after other conversions.
 */

#ifndef MP_FLAT_CVT_CONES_H
#define MP_FLAT_CVT_CONES_H

#include <cmath>
#include <limits>
#include <vector>
#include <cassert>

#include "mp/flat/constr_std.h"
#include "mp/flat/constr_keeper.h"
#include "mp/valcvt-link.h"


namespace mp {

/// Abstract converter for a single constraint.
template <class MCType, class Con>
class Convert1 { };


/** A class to store a ref to ModelConverter.
 *
 *  @param ModelConverter the converter type.
 */
template <class ModelConverter>
class MCKeeper {
public:
	/// Typedef MCType
	using MCType = ModelConverter;

	/// Constructor
	MCKeeper(MCType& mc) : mc_(mc) { }

protected:
	/// Retrieve the MC
	const MCType& MC() const { return mc_; }
	MCType& MC() { return mc_; }

private:
	MCType& mc_;
};


/**
 * Functor to convert various constraints to cones.
 *
 * @param MCType the main converter.
 */
template <class MCType>
class ConicConverter : public MCKeeper<MCType> {
public:
	virtual ~ConicConverter() { }

	/// Constructor
	ConicConverter(MCType& mc) :
		MCKeeper<MCType>(mc) { }

	/// Run conversions
	void Run() {
		if (MC().IfPassSOCPCones()) {
			Walk<QuadConRange>();
			Walk<QuadConLE>();
			Walk<QuadConGE>();

			Walk<LinConRange>();
			Walk<LinConLE>();
			Walk<LinConGE>();
    } else
      if (MC().IfPassQuadCon() &&
          (MC().GetNumberOfAddable((PowConstraint*)0)>0 ||
           MC().GetNumberOfAddable((AbsConstraint*)0)>0)) {
      // Still collect QCones expressed by 2-norms.
      // They are to be converted to quadratics.
      Walk<LinConRange>();
      Walk<LinConLE>();
      Walk<LinConGE>();
    }
	}


protected:
	/// Walk a single constraint type
	template <class Con>
	void Walk() {
		auto& ck = MC().GetConstraintKeeper((Con*)nullptr);
		ck.ForEachActive(Convert1<MCType, Con> { MC() });
	}

	/// Retrieve the MC
	using MCKeeper<MCType>::MC;
};


/** Converter for a single linear / quadratic
 *  inequality constraint.
 *  Generic, receives only body and sense.
 */
template <class MCType>
class Convert1QC : public MCKeeper<MCType> {
public:
	/// Constructor
	Convert1QC(MCType& mc) : MCKeeper<MCType>(mc) { }


protected:
  /// Characteristics of QP terms
  struct QPTermsTraits {
    int nDiffVars=0;
    int iDiffVars=-1;
    int x1=-1, x2=-1;
    double coef12=0.0;               // coefficient for x1*x2
    int nSamePos=0, nSameNeg=0;
    int iSamePos=-1, iSameNeg=-1;
  };

	/// DoRun. Body: quadratic.
	///
	/// Accept non-(+-1) coefficients.
	///
	/// @param body: quadratic constraint body.
	/// @param sens: -1 for <=, 1 for >=.
	/// @param rhs: constraint's rhs.
	///
	/// @return true iff converted.
	bool DoRun(const QuadAndLinTerms& body,
						 int sens, double rhs) {
		assert((sens==1 || sens==-1) && "sens 1 or -1 only");
    if (body.GetLinTerms().size()<=1) {
      QPTermsTraits qptt;
      if (Characterize(body.GetQPTerms(), qptt))
        return ProceedQCWithTraits(body, sens, rhs, qptt);
		}
		return false;
	}

  /// Characterize QP terms
  /// @return false iff not suitable
  bool Characterize(const QuadTerms& qpterms, QPTermsTraits& qptt) {
    for (auto i=qpterms.size(); i--;) {
      if (qpterms.var1(i) != qpterms.var2(i)) {
        if (qptt.nDiffVars++)   // Rotated cone: 2*x1*x2 >= ||xk||^2
          return false;
        qptt.x1 = qpterms.var1(i);
        qptt.x2 = qpterms.var2(i);
        if (MC().lb(qptt.x1)<0.0 || MC().lb(qptt.x2)<0.0)
          return false;              // Need nonnegative x1, x2
        qptt.coef12 = qpterms.coef(i);
        qptt.iDiffVars=i;
      } else {
        if (qpterms.coef(i)>0.0) {
          ++qptt.nSamePos;
          qptt.iSamePos=i;
        } else {
          ++qptt.nSameNeg;
          qptt.iSameNeg=i;
        }
      }
    }
    return true;
  }

  /// Proceed with a quadratic constraint provided traits.
  ///
  /// Cases of StdSOC -- in the "hacked form" x1^2 >= ||xk||^2:
  /// 1. x^2 >= ||y||^2
  /// 2. c^2 >= ||y||^2
  /// 3. x^2 >= ||y||^2 + c^2
  /// 4. 0   >= ||y||^2 + c^2  -- fail?
  ///
  /// Cases of RotSOC:
  /// 1. 2xy >= ||z||^2 [ + c^2 ]
  /// 2. 2 c1 y >= ||z||^2 [ + c2^2 ]
  ///
  /// @return true iff cone added.
  bool ProceedQCWithTraits(const QuadAndLinTerms& body,
                           int sens, double rhs,
                           const QPTermsTraits& qptt) {
    const auto& lint = body.GetLinTerms();
    assert(lint.size() <= 1);
    assert(qptt.nDiffVars <= 1);
    const auto& qpterms = body.GetQPTerms();
    if (qptt.nDiffVars || lint.size()) {       // Rotated cone?
      const auto fLT1 = lint.size()==1;        // Constant in the lhs?
      const auto fND1 = qptt.nDiffVars==1;
      if (fLT1 + fND1 == 1 &&
          rhs * sens >= 0.0) {                 // ... >= ... + rhs
        if (fLT1 &&                            // Check lin term
            (MC().lb( lint.var(0) ) < 0.0 ||
             lint.coef(0)*sens < 0.0))
          return false;
        if (fND1 && qptt.coef12*sens<0.0)
          return false;
        if ((sens>0 &&
             qptt.nSamePos==0 && qptt.nSameNeg)
            ||
            (sens<0 &&
             qptt.nSameNeg==0 && qptt.nSamePos))
          return AddRotatedQC(qpterms, lint, rhs,
                              qptt.iDiffVars);
      }
    } else {                          // Standard cone?
      if (0.0 <= rhs*sens &&          // It's x^2 >= ||y^2|| + |rhs|
          sens>0 && qptt.nSamePos==1 &&
          MC().lb(qpterms.var1(qptt.iSamePos))>=0.0)
        return AddStandardQC(qpterms, qptt.iSamePos, 0.0, rhs);
      if (0.0 <= rhs*sens &&
          sens<0 && qptt.nSameNeg==1 &&
          MC().lb(qpterms.var1(qptt.iSameNeg))>=0.0)
        return AddStandardQC(qpterms, qptt.iSameNeg, 0.0, rhs);
      if (0.0 >= rhs*sens &&          // |rhs| >= ||y^2||
          sens>0 && qptt.nSamePos==0 && qptt.nSameNeg)
        return AddStandardQC(qpterms, -1, rhs);
      if (0.0 >= rhs*sens &&
          sens<0 && qptt.nSameNeg==0 && qptt.nSamePos)
        return AddStandardQC(qpterms, -1, rhs);
    }
    return false;
  }

  /// DoRun.
  /// Body: linear, i.e., var1 (or const) >= var2 (or const).
  /// In particular, the following cases.
  ///
  /// Standard SOC.
  /// 1. x >= ||y [ | const ]||  (sqrt, abs)
  /// 2. const >= ||y [ | const ]||.
  ///
  /// Rotated SOC.
  /// 1. sqrt(k*x1*x2) >= ||z [ | const ]||.
  /// 2. sqrt(k * x2) >= ... (one var is const in lhs)
	///
  /// Considering const>=0.
	/// Accept non-(+-1) coefficients.
	///
	/// @param body: linear constraint body.
	/// @param sens: -1 for <=, 1 for >=.
	/// @param rhs: constraint's rhs.
	///
	/// @return true iff converted.
	bool DoRun(const LinTerms& lint,
						 int sens, double rhs) {
		assert((sens==1 || sens==-1) && "sens 1 or -1 only");
		if (2 == lint.size() &&                  // x>=y
				0.0 == std::fabs(rhs) &&             // rhs=0
				0.0 > lint.coef(0)*lint.coef(1)) {   // different signs
			// For sens>0 we have x>=y and iX=(index of pos coef)
			int iX = (lint.coef(1)>0.0);
			if (sens<0)
				iX = 1-iX;
			int iY = 1-iX;
			// Check if var(iY) := || vector-or-var ||
			if (auto rhs_args = CheckNorm2(lint.var(iY))) {
				if (auto lhs_args = CheckSqrtXnXmNonneg(lint.var(iX)))
					return ContinueRotatedSOC(lint, iX, iY,
																		lhs_args, rhs_args);
        return ContinueStdSOC(lint.coef(iX), lint.var(iX),
                              lint.coef(iY), rhs_args);
			}
    } else if (1 == lint.size() &&           // const>=y.
               0.0 >= rhs*sens &&            // either ... <= rhs
                                             // or ... >= rhs(<0)
               0.0 >= lint.coef(0)*sens) {   // In the other case,
                                             // it's -k*sqrt(...) <= rhs(>0), k>0,
                                             // which is always true TODO
      if (auto rhs_args = CheckNorm2(lint.var(0))) {
        return ContinueStdSOC(std::fabs(rhs),
                              MC().MakeFixedVar( 1.0 ),
                              lint.coef(0), rhs_args);
      }
    }
		return false;
	}

	/// Typedef for subexpression checkup result,
	/// whether it represents some part of an SOCP cone.
	struct ConeArgs {
		std::vector<double> coefs_;
		std::vector<int> vars_;
    double const_term = 0.0;
    double coef_extra = 0.0;
		/// Result vars of expressions to un-use.
		std::vector<int> res_vars_to_delete_;
		/// operator bool
		operator bool() const
		{ assert(check()); return !coefs_.empty(); }
		/// size()
		size_t size() const
		{ assert(check()); return coefs_.size(); }
		/// check()
		bool check() const { return coefs_.size()==vars_.size(); }
	};

	/// Check if the variable is defined by an expression
	/// representing ||.||.
	/// @return pair of (coef, var) vectors
	/// so that the cone is ... >= sqrt(sum{ (coef_i * var*i)^2 })
	ConeArgs CheckNorm2(int res_var) {
		if (MC().HasInitExpression(res_var)) {
			auto init_expr = MC().GetInitExpression(res_var);
			if (MC().template IsConInfoType<PowConstraint>(init_expr))
				return CheckNorm2_Pow(init_expr, res_var);
			if (MC().template IsConInfoType<AbsConstraint>(init_expr))
				return CheckNorm2_Abs(init_expr, res_var);
		}
		return {};
	}

	/// Check if the variable is defined by an expression
  /// representing sqrt(c1 * x1^2 + ... [ + const ]).
	template <class ConInfo>
	ConeArgs CheckNorm2_Pow(const ConInfo& ci, int res_var) {
		const auto& con_pow = MC().template
				GetConstraint<PowConstraint>(ci);
		const auto arg_pow = con_pow.GetArguments()[0];
		if (0.5 == con_pow.GetParameters()[0] &&     // sqrt(arg_pow)
				MC().HasInitExpression(arg_pow)) {
			auto ie_pow = MC().GetInitExpression(arg_pow);
			if (MC().template                          // arg_pow := QFC(...)
					IsConInfoType<QuadraticFunctionalConstraint>(ie_pow)) {
				const auto& con_qdc = MC().template
						GetConstraint<QuadraticFunctionalConstraint>(ie_pow);
				const auto& args_qdc = con_qdc.GetArguments();
        // Allowing a nonneg constant term by adding
				// an auxiliary variable^2
        if (0.0 <= args_qdc.constant_term() &&
						args_qdc.GetBody().GetLinTerms().empty()) {
					const auto& qpterms = args_qdc.GetBody().GetQPTerms();
					for (auto i=qpterms.size(); i--; ) {
						if (qpterms.coef(i) < 0.0 ||
								qpterms.var1(i) != qpterms.var2(i))
							return {};
					}
					ConeArgs result;
					result.coefs_ = qpterms.coefs();
					for (auto& c: result.coefs_)
						c = std::sqrt(c);
					result.vars_ = qpterms.vars1();
          result.const_term = args_qdc.constant_term();
					result.res_vars_to_delete_ = { res_var, arg_pow };
					return result;
				}
			}
		}
		return {};
	}

	/// Check if the variable is defined by an expression
	/// representing abs( c1*x1 ).
	template <class ConInfo>
	ConeArgs CheckNorm2_Abs(const ConInfo& ci, int res_var) {
		ConeArgs result;
		const auto& con_abs = MC().template
				GetConstraint<AbsConstraint>(ci);
		const auto arg_abs = con_abs.GetArguments()[0];
		result.coefs_ = { 1.0 };
		result.vars_ = { arg_abs };
		result.res_vars_to_delete_ = { res_var };
		return result;
	}

	/// Check if the variable is defined by sqrt(xN*xM) with
	/// xN, xM >= 0.
	ConeArgs CheckSqrtXnXmNonneg(int res_var) {
		ConeArgs result;
		if (auto pConPow = MC().template
				GetInitExpressionOfType<PowConstraint>(res_var)) {
			if (0.5 == pConPow->GetParameters()[0]) {     // sqrt(arg_pow)
				const auto arg_pow = pConPow->GetArguments()[0];
				if (auto pConQfc = MC().template
            GetInitExpressionOfType<           // arg_pow := QuadExpr
						QuadraticFunctionalConstraint>(arg_pow)) {
					const auto& qe = pConQfc->GetArguments();
					if (0.0 == std::fabs(qe.constant_term()) &&
							qe.GetBody().GetLinTerms().empty() &&
							1 == qe.GetBody().GetQPTerms().size()) {
						const auto& qpterms = qe.GetBody().GetQPTerms();
            if (qpterms.coef(0) >= 0.0 &&  // don't check if diff vars
								MC().lb(qpterms.var1(0)) >= 0.0 &&
								MC().lb(qpterms.var2(0)) >= 0.0) {
							result.coefs_ = {qpterms.coef(0), 1.0};
							result.vars_ = {qpterms.var1(0), qpterms.var2(0)};
							result.res_vars_to_delete_ = {res_var, arg_pow};
							return result;
						}
					}
        } else {        // It's just some variable
          if (MC().lb(arg_pow) >= 0.0) {
            result.coefs_ = {1.0};
            result.vars_ = {arg_pow};
            result.coef_extra = 1.0;   // for a new fixed var
            result.res_vars_to_delete_ = {res_var};
            return result;
          }
        }
			}
		}
		return result;
	}

	/// Continue processing a linear constraint x>=y,
  /// if y := ||.|| and x is sqrt(xN*xM) or sqrt(k*x).
	/// Rotated Qcone.
	/// @return true iff converted.
	bool ContinueRotatedSOC(
			const LinTerms& lint, int iX, int iY,
			const ConeArgs& lhs_args, const ConeArgs& rhs_args) {
    assert(2>=lhs_args.size());
    assert(1<=lhs_args.size());
    assert(rhs_args);
    const size_t constNZ =
        (0.0 != std::fabs(rhs_args.const_term));
    std::vector<double> c(2+rhs_args.size()+constNZ);
    std::vector<int> x(c.size());
		auto coefX_abs = std::fabs(lint.coef(iX));
		c[0] = 0.5 * lhs_args.coefs_[0] * coefX_abs;
		x[0] = lhs_args.vars_[0];
    if (2==lhs_args.size()) {   // x is sqrt(xN*xM)
      c[1] = lhs_args.coefs_[1] * coefX_abs;
      x[1] = lhs_args.vars_[1];
    } else {                    // x is sqrt(k*x)
      c[1] = lhs_args.coef_extra * coefX_abs;
      x[1] = MC().MakeFixedVar(1.0);
    }
		auto coefY_abs = std::fabs(lint.coef(iY));
		for (size_t iPush=0; iPush<rhs_args.size(); ++iPush) {
			x[iPush+2] = rhs_args.vars_[iPush];
			c[iPush+2] = coefY_abs * rhs_args.coefs_[iPush];
		}
    if (constNZ) {
      assert(0.0 < rhs_args.const_term);
      x.back() = MC().MakeFixedVar(1.0);
      c.back() = coefY_abs * std::sqrt(rhs_args.const_term);
    }
    for (auto r: lhs_args.res_vars_to_delete_)
			MC().DecrementVarUsage(r);
		for (auto r: rhs_args.res_vars_to_delete_)
			MC().DecrementVarUsage(r);
		MC().AddConstraint(
					RotatedQuadraticConeConstraint(
						std::move(x), std::move(c)));
		return true;
	}

  /// Continue processing a linear constraint cX*x >= cY*y,
	/// if y := ||.|| and x is not sqrt(xN*xM).
	/// Standard Qcone.
	/// @return true iff converted.
	bool ContinueStdSOC(
      double cX, int vX, double cY,
			const ConeArgs& rhs_args) {
    const size_t constNZ = (0.0 != std::fabs(rhs_args.const_term));
    std::vector<int> x(rhs_args.size()+1+constNZ);
    std::vector<double> c(rhs_args.size()+1+constNZ);
    x[0] = vX;
    c[0] = std::fabs(cX);
    auto coefY_abs = std::fabs(cY);
		for (size_t iPush=0; iPush<rhs_args.size(); ++iPush) {
			x[iPush+1] = rhs_args.vars_[iPush];
			c[iPush+1] = coefY_abs * rhs_args.coefs_[iPush];
		}
    if (constNZ) {
      assert(0.0 < rhs_args.const_term);
      x.back() = MC().MakeFixedVar(1.0);
      c.back() = coefY_abs * std::sqrt(rhs_args.const_term);
    }
		for (auto r: rhs_args.res_vars_to_delete_)
			MC().DecrementVarUsage(r);
		MC().AddConstraint(
					QuadraticConeConstraint(
						std::move(x), std::move(c)));
		return true;
	}


	/// Add rotated cone from pure-quadratic constraint
  bool AddRotatedQC(const QuadTerms& qpterms,
                    const LinTerms& lint,
                    double rhs, int iDiffVars) {
    assert(lint.size() + (iDiffVars>=0) == 1);
    const auto rhsNon0 = (std::fabs(rhs)!=0.0);
    std::vector<int> x(qpterms.size()+1+lint.size()+rhsNon0);
    std::vector<double> c(qpterms.size()+1+lint.size()+rhsNon0);
		size_t iPush=1;
    if (lint.size()) {
      x[0] = lint.var(0);
      x[1] = MC().MakeFixedVar( 1.0 );
      c[0] = std::fabs(lint.coef(0)) / 2.0;
      c[1] = 1.0;      // it's 2xy >= z^2 + ...
    }
		for (int i=0; i<(int)qpterms.size(); ++i) {
			if (i==iDiffVars) {
				x[0] = qpterms.var1(i);
				x[1] = qpterms.var2(i);
				c[0] = std::fabs(qpterms.coef(i));
        c[1] = 0.5;           // it's 2xy >= z^2 + ...
			} else {
				++iPush;
				x.at(iPush) = qpterms.var1(i);
				c.at(iPush) = std::sqrt(std::fabs(qpterms.coef(i)));
			}
		}
    if (rhsNon0) {            // Add  ... >= ... + |rhs|
      ++iPush;
      x.at(iPush) = MC().MakeFixedVar( 1.0 );
      c.at(iPush) = std::sqrt( std::fabs(rhs) );
    }
    assert(iPush == x.size()-1);
		MC().AddConstraint(
					RotatedQuadraticConeConstraint(
						std::move(x), std::move(c)));
		return true;
	}

  /// Add standard cone from pure-quadratic constraint.
  /// @param iSame1: if <0, use new fixed variable.
  /// @param rhs: square of the new value if iSame1<0,
  ///         i.e., the constraint is rhs >= ....
  /// @param rhs2: ... >= ... + rhs2
  bool AddStandardQC(const QuadTerms& qpterms,
                     int iSame1,
                     double rhs=0.0, double rhs2=0.0) {
    const size_t fRhs2 = (0.0!=std::fabs(rhs2));
    const size_t fNewVar = (iSame1<0) || fRhs2;
    std::vector<int> x(qpterms.size() + fNewVar);
    std::vector<double> c(qpterms.size() + fNewVar);
    if (fNewVar && !fRhs2) {
      c[0] = std::sqrt(std::fabs(rhs));
      x[0] = MC().MakeFixedVar( 1.0 );
    }
    size_t iPush=0;
		for (int i=0; i<(int)qpterms.size(); ++i) {
			if (i==iSame1) {
				x[0] = qpterms.var1(i);
				c[0] = std::sqrt(std::fabs(qpterms.coef(i)));
			} else {
				++iPush;
				x.at(iPush) = qpterms.var1(i);
				c.at(iPush) = std::sqrt(std::fabs(qpterms.coef(i)));
			}
		}
    if (fNewVar && fRhs2) {
      ++iPush;
      c.at(iPush) = std::sqrt(std::fabs(rhs2));
      x.at(iPush) = MC().MakeFixedVar( 1.0 );
    }
    assert(iPush == c.size()-1);
    MC().AddConstraint(
					QuadraticConeConstraint(
						std::move(x), std::move(c)));
		return true;
	}

	/// Retrieve the MC
	using MCKeeper<MCType>::MC;
};


/** Converter for a single algebraic inequality constraint.
 *  Receives range constraints and proceeds
 *  if they are inequalities.
 */
template <class MCType, class Body>
class Convert1<MCType,
		AlgebraicConstraint<Body, AlgConRange> > :
		public Convert1QC<MCType> {
public:
	/// Constructor
	Convert1(MCType& mc) : Convert1QC<MCType>(mc) { }

	/// Constraint type
	using ConType = AlgebraicConstraint<Body, AlgConRange>;

	/// Run
	bool operator()(const ConType& con, int i) {
		bool lbf = con.lb() >= MC().PracticallyMinusInf();
		bool ubf = con.ub() <= MC().PracticallyInf();
		if (lbf+ubf == 1) {
			// Autolink is to pass back duals... Simplify?
			// E.g., automatically in ForEachActive?
			pre::AutoLinkScope<MCType> auto_link_scope{
				MC(), MC().template SelectValueNode<ConType>(i)
			};
			return
					Convert1QC<MCType>::DoRun(con.GetBody(),
																lbf ? 1 : -1,
																lbf ? con.lb() : con.ub());
		}
		return false;
	}

	/// Retrieve the MC
	using MCKeeper<MCType>::MC;
};


/** Converter for a single algebraic inequality constraint.
*   Proper inequality constraints.
*/
template <class MCType, class Body, int sens>
class Convert1<MCType,
		AlgebraicConstraint< Body, AlgConRhs<sens> > > :
		public Convert1QC<MCType> {
public:
	/// Constructor
	Convert1(MCType& mc) : Convert1QC<MCType>(mc) { }

	/// The constraint type
	using ConType =
		AlgebraicConstraint< Body, AlgConRhs<sens> >;

	/// Run
	bool operator()(const ConType& con, int i) {
		static_assert (sens==1 || sens==-1, "sens 1 or -1 only");
		pre::AutoLinkScope<MCType> auto_link_scope{
			MC(), MC().template SelectValueNode<ConType>(i)
		};
		return
				Convert1QC<MCType>::DoRun(con.GetBody(),
																	sens, con.rhs());
	}

	/// Retrieve the MC
	using MCKeeper<MCType>::MC;
};

}  // namespace mp

#endif // MP_FLAT_CVT_CONES_H
