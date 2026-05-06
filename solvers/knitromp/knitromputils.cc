#include <functional>

#include "mp/format.h"
#include "knitromputils.h"

namespace mp {

    // Recursive function to print expression tree to MemoryWriter
    void ExpressionData::FormatExpression(fmt::MemoryWriter& w, const ExpressionData& expr,
        std::function<void(int, fmt::MemoryWriter& w)> getVarName) {

        auto opName = [](ExpressionData::OpType op) -> const char* {
            switch (op) {
            case ExpressionData::VAR: return "VAR";
            case ExpressionData::CONST: return "CONST";
            case ExpressionData::ADD: return "+";
            case ExpressionData::SUB: return "-";
            case ExpressionData::MUL: return "*";
            case ExpressionData::DIV: return "/";
            case ExpressionData::SIN: return "sin";
            case ExpressionData::COS: return "cos";
            case ExpressionData::TAN: return "tan";
            case ExpressionData::LOG: return "log";
            case ExpressionData::EXP: return "exp";
            case ExpressionData::POW: return "^";
            case ExpressionData::SQRT: return "sqrt";
            case ExpressionData::ASIN: return "asin";
            case ExpressionData::ACOS: return "acos";
            case ExpressionData::ATAN: return "atan";
            case ExpressionData::SINH: return "sinh";
            case ExpressionData::COSH: return "cosh";
            case ExpressionData::TANH: return "tanh";
            case ExpressionData::ASINH: return "asinh";
            case ExpressionData::ACOSH: return "acosh";
            case ExpressionData::ATANH: return "atanh";
            case ExpressionData::ABS: return "abs";
            default: return "UNKNOWN";
            }
            };

        std::function<void(const ExpressionData&)> formatExpr =
            [&](const ExpressionData& e) -> void {

            switch (e.op_) {
            case ExpressionData::VAR:
                getVarName(e.varIndex_, w);
                break;

            case ExpressionData::CONST:
                w << e.constValue_;
                break;


            case ExpressionData::MUL:
            case ExpressionData::ADD:
            case ExpressionData::SUB:
            case ExpressionData::DIV:
            case ExpressionData::POW:
                w << "(";
                if (e.left_) {
                    formatExpr(*e.left_);
                }
                else {
                    w << "?";
                }
                w << " " << opName(e.op_) << " ";
                if (e.right_) {
                    formatExpr(*e.right_);
                }
                else {
                    w << "?";
                }
                w << ")";
                break;

            default: {  // Unary operations
                w << opName(e.op_) << "(";
                if (e.left_) {
                    formatExpr(*e.left_);
                }
                else {
                    w << "?";
                }
                w << ")";
                break;
            }
            }
            };

        formatExpr(expr);
        if (expr.linearCoeffs().size() > 0) {
            w << " ";
            for (size_t i = 0; i < expr.linearCoeffs().size(); i++) {
                PrintConstraint::PrintCoef(expr.linearCoeffs()[i], w, false);
                getVarName(expr.linearIndices()[i], w);
            }
        }
    }



}