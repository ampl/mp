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


/** Converter for a single quadratic inequality constraint.
 *  Generic, receives only body (assumes rhs=0) and sense.
 */
template <class MCType>
class Convert1QC : public MCKeeper<MCType> {
public:
	/// Constructor
	Convert1QC(MCType& mc) : MCKeeper<MCType>(mc) { }


protected:
	/// DoRun. Body: quadratic.
	///
	/// Currently considering rhs=0.
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
		if (body.GetLinTerms().empty() &&
				0.0 == std::fabs(rhs)) {       // rhs=0
			const auto& qpterms = body.GetQPTerms();
			int nDiffVars=0;
			int iDiffVars=-1;
			int x1=-1, x2=-1;
			double coef12=0.0;               // coefficient for x1*x2
			int nSamePos=0, nSameNeg=0;
			int iSamePos=-1, iSameNeg=-1;
			for (auto i=qpterms.size(); i--;) {
				if (qpterms.var1(i) != qpterms.var2(i)) {
					if (nDiffVars++)   // Rotated cone: 2*x1*x2 >= ||xk||^2
						return false;
					x1 = qpterms.var1(i);
					x2 = qpterms.var2(i);
					if (MC().lb(x1)<0.0 || MC().lb(x2)<0.0)
						return false;              // Need nonnegative x1, x2
					coef12 = qpterms.coef(i);
					iDiffVars=i;
				} else {
					if (qpterms.coef(i)>0.0) {
						++nSamePos;
						iSamePos=i;
					} else {
						++nSameNeg;
						iSameNeg=i;
					}
				}
			}
			if (nDiffVars) {                 // Rotated cone?
				if ((sens>0 && coef12>0.0 && nSamePos==0 && nSameNeg)
						||
						(sens<0 && coef12<0.0 && nSameNeg==0 && nSamePos))
					return AddRotatedQC(qpterms, iDiffVars);
			} else {                         // Standard cone?
				if (sens>0 && nSamePos==1 && nSameNeg &&
						MC().lb(qpterms.var1(iSamePos))>=0.0)
												// In the "hacked form" x1^2 >= ||xk||^2
					return AddStandardQC(qpterms, iSamePos);
				if (sens<0 && nSameNeg==1 && nSamePos &&
						MC().lb(qpterms.var1(iSameNeg))>=0.0)
					return AddStandardQC(qpterms, iSameNeg);
			}
		}
		return false;
	}

	/// DoRun. Body: linear.
	///
	/// Currently considering rhs=0.
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
			// Check if iY = || vector-or-var ||
			if (auto rhs_args = CheckNorm2(lint.var(iY))) {
				std::vector<int> x(rhs_args.size()+1);
				std::vector<double> c(rhs_args.size()+1);
				x[0] = lint.var(iX);
				c[0] = std::fabs(lint.coef(iX));
				auto coefY_abs = std::fabs(lint.coef(iY));
				for (size_t iPush=0; iPush<rhs_args.size(); ++iPush) {
					x[iPush+1] = rhs_args.vars_[iPush];
					c[iPush+1] = coefY_abs * rhs_args.coefs_[iPush];
				}
				// Un-use y together with any subexpression result vars.
				rhs_args.res_vars_to_delete_.push_back(lint.var(iY));
				for (auto r: rhs_args.res_vars_to_delete_)
					MC().DecrementVarUsage(r);
				MC().AddConstraint(
							QuadraticConeConstraint(
								std::move(x), std::move(c)));
			}
		}
		return false;
	}

	/// Typedef for subexpression checkup result,
	/// whether it represents the RHS of an SOCP cone.
	/// If yes, populates members so that the RHS is
	/// sqrt( sum { coef_[i]*vars_[i] } )
	struct ConeRHSArgs {
		std::vector<double> coefs_;
		std::vector<int> vars_;
		/// Result vars of expressions to un-use.
		std::vector<int> res_vars_to_delete_;
		/// operator bool
		operator bool() const { return ! coefs_.empty(); }
		/// size()
		size_t size() const { return coefs_.size(); }
	};

	/// Check if the variable is defined by an expression
	/// representing ||.||.
	/// @return pair of (coef, var) vectors
	/// so that the cone is ... >= sqrt(sum{ (coef_i * var*i)^2 })
	ConeRHSArgs	CheckNorm2(int res_var) {
		if (MC().HasInitExpression(res_var)) {
			auto init_expr = MC().GetInitExpression(res_var);
			if (MC().template IsConInfoType<PowConstraint>(init_expr))
				return CheckNorm2_Pow(init_expr);
			if (MC().template IsConInfoType<AbsConstraint>(init_expr))
				return CheckNorm2_Abs(init_expr);
		}
		return {};
	}

	/// Check if the variable is defined by an expression
	/// representing sqrt(c1 * x1^2 + ...).
	template <class ConInfo>
	ConeRHSArgs	CheckNorm2_Pow(const ConInfo& ci) {
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
				// We could allow a nonneg constant term by adding
				// an auxiliary variable^2
				if (0.0 == std::fabs(args_qdc.constant_term()) &&
						args_qdc.GetBody().GetLinTerms().empty()) {
					const auto& qpterms = args_qdc.GetBody().GetQPTerms();
					for (auto i=qpterms.size(); i--; ) {
						if (qpterms.coef(i) < 0.0 ||
								qpterms.var1(i) != qpterms.var2(i))
							return {};
					}
					ConeRHSArgs result;
					result.coefs_ = qpterms.coefs();
					for (auto& c: result.coefs_)
						c = std::sqrt(c);
					result.vars_ = qpterms.vars1();
					result.res_vars_to_delete_ = { arg_pow };
					return result;
				}
			}
		}
		return {};
	}

	/// Check if the variable is defined by an expression
	/// representing abs( c1*x1 ).
	template <class ConInfo>
	ConeRHSArgs	CheckNorm2_Abs(const ConInfo& ) {
		ConeRHSArgs result;
		// Probably better linearize Abs()
		return result;
	}

	/// Add rotated cone
	bool AddRotatedQC(const QuadTerms& qpterms, int iDiffVars) {
		std::vector<int> x(qpterms.size()+1);
		std::vector<double> c(qpterms.size()+1);
		size_t iPush=1;
		for (int i=0; i<(int)qpterms.size(); ++i) {
			if (i==iDiffVars) {
				x[0] = qpterms.var1(i);
				x[1] = qpterms.var2(i);
				c[0] = std::fabs(qpterms.coef(i));
				c[1] = 0.5;      // it's 2xy >= z^2 + ...
			} else {
				++iPush;
				x.at(iPush) = qpterms.var1(i);
				c.at(iPush) = std::sqrt(std::fabs(qpterms.coef(i)));
			}
		}
		MC().AddConstraint(
					RotatedQuadraticConeConstraint(
						std::move(x), std::move(c)));
		return true;
	}

	/// Add standard cone
	bool AddStandardQC(const QuadTerms& qpterms, int iSame1) {
		std::vector<int> x(qpterms.size());
		std::vector<double> c(qpterms.size());
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
