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
		if (MC().IfPassQuadCones()) {
			Walk<QuadConRange>();
			Walk<QuadConLE>();
			Walk<QuadConGE>();
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
	/// DoRun.
	///
	/// Currently considering rhs=0.
	/// Might do non-(+-1) coefficients.
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

	/// Add rotated cone
	bool AddRotatedQC(const QuadTerms& qpterms, int iDiffVars) {
		std::vector<int> x(qpterms.size()+1);
		std::vector<double> c;
		c.reserve(qpterms.size()-1);
		size_t iPush=1;
		for (auto i=(int)qpterms.size(); i--; ) {
			if (i==iDiffVars) {
				x[0] = qpterms.var1(i);
				x[1] = qpterms.var2(i);
			} else {
				x.at(++iPush) = qpterms.var1(i);
				c.push_back(std::fabs(qpterms.coef(i)));
			}
		}
		MC().AddConstraint(RotatedQuadraticConeConstraint(x, c));
		return true;
	}

	/// Add standard cone
	bool AddStandardQC(const QuadTerms& qpterms, int iSame1) {
		std::vector<int> x(qpterms.size());
		std::vector<double> c;
		c.reserve(qpterms.size()-1);
		size_t iPush=0;
		for (auto i=(int)qpterms.size(); i--; ) {
			if (i==iSame1) {
				x[0] = qpterms.var1(i);
			} else {
				x.at(++iPush) = qpterms.var1(i);
				c.push_back(std::fabs(qpterms.coef(i)));
			}
		}
		MC().AddConstraint(QuadraticConeConstraint(x, c));
		return true;
	}

	/// Retrieve the MC
	using MCKeeper<MCType>::MC;
};


/** Converter for a single quadratic inequality constraint.
 *  Receives range constraints and proceeds
 *  if they are inequalities.
 */
template <class MCType>
class Convert1<MCType, QuadConRange> : public Convert1QC<MCType> {
public:
	/// Constructor
	Convert1(MCType& mc) : Convert1QC<MCType>(mc) { }

	/// Run
	bool operator()(const QuadConRange& con, int i) {
		bool lbf = con.lb() >= MC().PracticallyMinusInf();
		bool ubf = con.ub() <= MC().PracticallyInf();
		if (lbf+ubf == 1) {
			// Autolink is to pass back duals... Simplify?
			// E.g., automatically in ForEachActive?
			auto& ck = MC().GetConstraintKeeper((QuadConRange*)nullptr);
			pre::AutoLinkScope<MCType> auto_link_scope{
				MC(), ck.SelectValueNodeRange(i)
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


/** Converter for a single quadratic inequality constraint.
*   Proper inequality constraints.
*/
template <class MCType, int sens>
class Convert1<MCType, QuadConRhs<sens> > :
		public Convert1QC<MCType> {
public:
	/// Constructor
	Convert1(MCType& mc) : Convert1QC<MCType>(mc) { }

	/// The constraint type
	using ConType = QuadConRhs<sens>;

	/// Run
	bool operator()(const ConType& con, int i) {
		static_assert (sens==1 || sens==-1, "sens 1 or -1 only");
		auto& ck = MC().GetConstraintKeeper((ConType*)nullptr);
		pre::AutoLinkScope<MCType> auto_link_scope{
			MC(), ck.SelectValueNodeRange(i)
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
