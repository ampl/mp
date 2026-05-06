#ifndef KNITROMPUTILS_H
#define KNITROMPUTILS_H


#include <vector>
#include <set>

#include "mp/format.h"
#include "knitro.h" // for KN_INFINITY

namespace mp {

    class PrintConstraint {
        double lhs, rhs;
        std::vector<int> linearIndices;
        std::vector<double> linearCoeffs;
        std::vector<int> quadVars1;
        std::vector<int> quadVars2;
        std::vector<double> quadCoeffs;
    public:
        PrintConstraint() {}
		PrintConstraint(double lhs, double rhs,const std::vector<int> &linearIndices, const std::vector<double> &linearCoeffs)
			: lhs(lhs), rhs(rhs), linearIndices(linearIndices), linearCoeffs(linearCoeffs)
			 {}

        PrintConstraint(double lhs, double rhs, const std::vector<int> &linearIndices, const std::vector<double> &linearCoeffs,
            const std::vector<int> &quadVars1, const std::vector<int> &quadVars2, const std::vector<double> &quadCoeffs)
            : lhs(lhs), rhs(rhs), linearIndices(linearIndices), linearCoeffs(linearCoeffs),
            quadVars1(quadVars1), quadVars2(quadVars2), quadCoeffs(quadCoeffs) {}


        static void PrintCoef(double c, fmt::MemoryWriter &w, bool first)  {
            if (c > 0) {
                if (!first) w << " + ";
                // Omit '1*' for positive coefficient 1
                if (c != 1)
                    w << c << "*";
            }
            else { // Coefficient is negative
                // If the coefficient is -1, just print ' - ', otherwise print the coefficient
                if (c == -1) w << " - ";
                else  w << c << "*";
            }
		}
        void printCoeffs(fmt::MemoryWriter &w, std::function<void(int, fmt::MemoryWriter& w)> getVarName) const {
            if (linearIndices.size() > 0)
            {
                for (size_t i = 0; i < linearIndices.size(); i++) {
					PrintCoef(linearCoeffs[i], w, i == 0);
					getVarName(linearIndices[i], w);
                }
            }
            if (quadVars1.size() > 0) {
                for (size_t i = 0; i < quadVars1.size(); i++) {
                    PrintCoef(quadCoeffs[i], w, i == 0);
                    getVarName(quadVars1[i], w);
                    w << " * ";
                    getVarName(quadVars2[i], w);
                }
            }
           
        }

        void print(fmt::MemoryWriter& w, std::function<void(int, fmt::MemoryWriter& w)> getVarName) const {
            bool eq = rhs == lhs;

			if (!eq && (lhs > -KN_INFINITY)) {
                w << lhs << " <= ";
            }
            printCoeffs(w, getVarName);
            if (!eq && ( rhs< KN_INFINITY)) {
                w << " <= ";
                w << rhs;
            }
            if (eq)
            {
                w << " = ";
                w << rhs;
            }
            w << "\n";
        }

        void printObjective(bool maximize, fmt::MemoryWriter &w, std::function<void(int, fmt::MemoryWriter& w)> getVarName) const {
            printCoeffs(w, getVarName);
            w << "\n";
        }


    };


    class ExpressionData {
    public:
        enum OpType {
            VAR,      // Variable reference
            CONST,    // Constant value
            ADD, SUB, POW, MUL, DIV,  // Binary arithmetic
            SIN, COS, TAN, ASIN, ACOS, ATAN,  // Trig
            SINH, COSH, TANH, ASINH, ACOSH, ATANH,  // Hyperbolic
            LOG, EXP, SQRT, ABS  // Other functions

        };


       


        void reserveLinear(int size) {
            linearCoeffs_.reserve(size);
            linearIndices_.reserve(size);
            empty_ = false;
        }
        void addLinear(int index, double coeff) {
            linearCoeffs_.push_back(coeff);
            linearIndices_.push_back(index);
            empty_ = false;
        }


        static ExpressionData BinaryOp(OpType type,
            const ExpressionData& lhs, const ExpressionData& rhs)
        {
            assert((type >= ADD) && (type <= DIV));
            ExpressionData result;
            result.op_ = type;
            result.left_ = std::make_shared<ExpressionData>(lhs);
            result.right_ = std::make_shared<ExpressionData>(rhs);
            result.usedVars_ = lhs.usedVars_;
            result.usedVars_.insert(rhs.usedVars_.begin(), rhs.usedVars_.end());
            result.empty_ = false;
            return result;

        }
        // Expression operations
        static ExpressionData Add(const ExpressionData& lhs, const ExpressionData& rhs) {
            return BinaryOp(ExpressionData::ADD, lhs, rhs);
        }

        static ExpressionData Multiply(const ExpressionData& lhs, const ExpressionData& rhs) {
            return BinaryOp(ExpressionData::MUL, lhs, rhs);
        }

        // Unary operations
        static ExpressionData MakeUnaryExpr(ExpressionData::OpType op, const ExpressionData& arg) {
            ExpressionData result;
            result.op_ = op;
            result.left_ = std::make_shared<ExpressionData>(arg);
            result.usedVars_ = arg.usedVars_;
            result.empty_ = false;
            return result;
        }

        /// Make a constant expression.
        static ExpressionData MakeConstantExpr(double v) {
            ExpressionData p;
            p.op_ = ExpressionData::CONST;
            p.constValue_ = v;
            p.empty_ = false;
            return p;
        }

        /// Make an empty expression.
        static ExpressionData MakeEmptyExpr() {
            return ExpressionData();
        }

        /// Make an expression representing variable \a v.
        static ExpressionData MakeVarExpr(int v) {
            ExpressionData p;
            p.op_ = ExpressionData::VAR;
            p.varIndex_ = v;
            p.usedVars_.insert(v);
            p.empty_ = false;
            return p;
        }

        static void ExpressionData::FormatExpression(fmt::MemoryWriter& w, const ExpressionData& expr,
            std::function<void(int, fmt::MemoryWriter& w)> getVarName);
        std::vector<int> linearIndices() const { return linearIndices_; }
        std::vector<double> linearCoeffs() const { return linearCoeffs_; }
        bool isEmpty() const { return empty_; }

        struct HessInfo {
            std::set<int> gradVars;
            std::set<std::pair<int, int>> hessPairs; // lower triangle: first >= second
        };

        static HessInfo computeHessianSparsityInfo(const ExpressionData& e)
        {
            HessInfo info;

            switch (e.op_) {
            case ExpressionData::CONST:
                return info;

            case ExpressionData::VAR:
                info.gradVars.insert(e.varIndex_);
                return info;

            case ExpressionData::ADD:
            case ExpressionData::SUB: {
                auto L = computeHessianSparsityInfo(*e.left_);
                auto R = computeHessianSparsityInfo(*e.right_);

                info.gradVars = L.gradVars;
                info.gradVars.insert(R.gradVars.begin(), R.gradVars.end());

                info.hessPairs = L.hessPairs;
                info.hessPairs.insert(R.hessPairs.begin(), R.hessPairs.end());
                return info;
            }

            case ExpressionData::MUL: {
                auto L = computeHessianSparsityInfo(*e.left_);
                auto R = computeHessianSparsityInfo(*e.right_);

                info.gradVars = L.gradVars;
                info.gradVars.insert(R.gradVars.begin(), R.gradVars.end());

                info.hessPairs = L.hessPairs;
                info.hessPairs.insert(R.hessPairs.begin(), R.hessPairs.end());

                AddCrossPairs(info.hessPairs, L.gradVars, R.gradVars);
                return info;
            }

            case ExpressionData::DIV:
            case ExpressionData::POW: {
                auto L = computeHessianSparsityInfo(*e.left_);
                auto R = computeHessianSparsityInfo(*e.right_);

                info.gradVars = L.gradVars;
                info.gradVars.insert(R.gradVars.begin(), R.gradVars.end());

                info.hessPairs = L.hessPairs;
                info.hessPairs.insert(R.hessPairs.begin(), R.hessPairs.end());

                // Conservative: division/power can couple all argument variables
                AddAllPairs(info.hessPairs, info.gradVars);
                return info;
            }

            case ExpressionData::SIN:
            case ExpressionData::COS:
            case ExpressionData::TAN:
            case ExpressionData::ASIN:
            case ExpressionData::ACOS:
            case ExpressionData::ATAN:
            case ExpressionData::SINH:
            case ExpressionData::COSH:
            case ExpressionData::TANH:
            case ExpressionData::ASINH:
            case ExpressionData::ACOSH:
            case ExpressionData::ATANH:
            case ExpressionData::LOG:
            case ExpressionData::EXP:
            case ExpressionData::SQRT: {
                auto A = computeHessianSparsityInfo(*e.left_);

                info.gradVars = A.gradVars;
                info.hessPairs = A.hessPairs;

                // f(g(x)) has Hessian:
                // f''(g) grad(g) grad(g)^T + f'(g) Hessian(g)
                AddAllPairs(info.hessPairs, A.gradVars);
                return info;
            }

            default:
                throw std::runtime_error("Unknown expression op in Hessian sparsity");
            }
        }


		OpType op() const { return op_; }
		int varIndex() const { assert(op_ == VAR); return varIndex_; }
		double constValue() const { assert(op_ == CONST); return constValue_; }
        std::shared_ptr<ExpressionData> left() const {
            return left_;
        };
        std::shared_ptr<ExpressionData> right() const {
			return right_;
        }
        const std::set<int> &usedVars() const {
			return usedVars_;
        }
    private:
        std::vector<int> linearIndices_;
        std::vector<double> linearCoeffs_;
        OpType op_;
        double constValue_ = 0.0;
        int varIndex_ = -1;
        std::shared_ptr<ExpressionData> left_;
        std::shared_ptr<ExpressionData> right_;
        std::set<int> usedVars_;  // All variables used in this expression tree
        bool empty_ = true;
        static std::pair<int, int> MakeLowerPair(int a, int b) {
            if (a < b)
                std::swap(a, b);
            return { a, b };
        }

        static void AddAllPairs(
            std::set<std::pair<int, int>>& out,
            const std::set<int>& vars)
        {
            for (int a : vars)
                for (int b : vars)
                    out.insert(MakeLowerPair(a, b));
        }

        static void AddCrossPairs(
            std::set<std::pair<int, int>>& out,
            const std::set<int>& aVars,
            const std::set<int>& bVars)
        {
            for (int a : aVars)
                for (int b : bVars)
                    out.insert(MakeLowerPair(a, b));
        }
    };



   
} // namespace mp

#endif // KNITROMPUTILS_H
