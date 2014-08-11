/*
 Base types and functions related to the optimization problem.

 Copyright (C) 2014 AMPL Optimization Inc

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

#include "mp/problem-base.h"
#include "precedence.h"

namespace mp {

const expr::Info expr::Info::INFO[N_OPS] = {
    {expr::BINARY,           prec::ADDITIVE,       "+"},  // OPPLUS
    {expr::BINARY,           prec::ADDITIVE,       "-"},  // OPMINUS
    {expr::BINARY,           prec::MULTIPLICATIVE, "*"},  // OPMULT
    {expr::BINARY,           prec::MULTIPLICATIVE, "/"},  // OPDIV
    {expr::BINARY,           prec::MULTIPLICATIVE, "mod"},  // OPREM
    {expr::BINARY,           prec::EXPONENTIATION, "^"},  // OPPOW
    {expr::BINARY,           prec::ADDITIVE,       "less"},  // OPLESS
    {expr::UNKNOWN,          prec::UNKNOWN,        "unknown"},
    {expr::UNKNOWN,          prec::UNKNOWN,        "unknown"},
    {expr::UNKNOWN,          prec::UNKNOWN,        "unknown"},
    {expr::UNKNOWN,          prec::UNKNOWN,        "unknown"},
    {expr::VARARG,           prec::CALL,           "min"},  // MINLIST
    {expr::VARARG,           prec::CALL,           "max"},  // MAXLIST
    {expr::UNARY,            prec::CALL,           "floor"},  // FLOOR
    {expr::UNARY,            prec::CALL,           "ceil"},  // CEIL
    {expr::UNARY,            prec::CALL,           "abs"},  // ABS
    {expr::UNARY,            prec::UNARY,          "unary -"},  // OPUMINUS
    {expr::UNKNOWN,          prec::UNKNOWN,        "unknown"},
    {expr::UNKNOWN,          prec::UNKNOWN,        "unknown"},
    {expr::UNKNOWN,          prec::UNKNOWN,        "unknown"},
    {expr::BINARY_LOGICAL,   prec::LOGICAL_OR,     "||"},  // OPOR
    {expr::BINARY_LOGICAL,   prec::LOGICAL_AND,    "&&"},  // OPAND
    {expr::RELATIONAL,       prec::RELATIONAL,     "<"},  // LT
    {expr::RELATIONAL,       prec::RELATIONAL,     "<="},  // LE
    {expr::RELATIONAL,       prec::RELATIONAL,     "="},  // EQ
    {expr::UNKNOWN,          prec::UNKNOWN,        "unknown"},
    {expr::UNKNOWN,          prec::UNKNOWN,        "unknown"},
    {expr::UNKNOWN,          prec::UNKNOWN,        "unknown"},
    {expr::RELATIONAL,       prec::RELATIONAL,     ">="},  // GE
    {expr::RELATIONAL,       prec::RELATIONAL,     ">"},  // GT
    {expr::RELATIONAL,       prec::RELATIONAL,     "!="},  // NE
    {expr::UNKNOWN,          prec::UNKNOWN,        "unknown"},
    {expr::UNKNOWN,          prec::UNKNOWN,        "unknown"},
    {expr::UNKNOWN,          prec::UNKNOWN,        "unknown"},
    {expr::NOT,              prec::NOT,            "!"},  // OPNOT
    {expr::IF,               prec::CONDITIONAL,    "if"},  // OPIFnl
    {expr::UNKNOWN,          prec::UNKNOWN,        "unknown"},
    {expr::UNARY,            prec::CALL,           "tanh"},  // OP_tanh
    {expr::UNARY,            prec::CALL,           "tan"},  // OP_tan
    {expr::UNARY,            prec::CALL,           "sqrt"},  // OP_sqrt
    {expr::UNARY,            prec::CALL,           "sinh"},  // OP_sinh
    {expr::UNARY,            prec::CALL,           "sin"},  // OP_sin
    {expr::UNARY,            prec::CALL,           "log10"},  // OP_log10
    {expr::UNARY,            prec::CALL,           "log"},  // OP_log
    {expr::UNARY,            prec::CALL,           "exp"},  // OP_exp
    {expr::UNARY,            prec::CALL,           "cosh"},  // OP_cosh
    {expr::UNARY,            prec::CALL,           "cos"},  // OP_cos
    {expr::UNARY,            prec::CALL,           "atanh"},  // OP_atanh
    {expr::BINARY,           prec::CALL,           "atan2"},  // OP_atan2
    {expr::UNARY,            prec::CALL,           "atan"},  // OP_atan
    {expr::UNARY,            prec::CALL,           "asinh"},  // OP_asinh
    {expr::UNARY,            prec::CALL,           "asin"},  // OP_asin
    {expr::UNARY,            prec::CALL,           "acosh"},  // OP_acosh
    {expr::UNARY,            prec::CALL,           "acos"},  // OP_acos
    {expr::SUM,              prec::ITERATIVE,      "sum"},  // OPSUMLIST
    {expr::BINARY,           prec::MULTIPLICATIVE, "div"},  // OPintDIV
    {expr::BINARY,           prec::CALL,           "precision"},  // OPprecision
    {expr::BINARY,           prec::CALL,           "round"},  // OPround
    {expr::BINARY,           prec::CALL,           "trunc"},  // OPtrunc
    {expr::COUNT,            prec::CALL,           "count"},  // OPCOUNT
    {expr::NUMBEROF,         prec::CALL,           "numberof"},  // OPNUMBEROF
    // OPNUMBEROFs - not supported yet
    {expr::UNKNOWN,          prec::UNKNOWN,        "string numberof"},
    {expr::LOGICAL_COUNT,    prec::CALL,           "atleast"},  // OPATLEAST
    {expr::LOGICAL_COUNT,    prec::CALL,           "atmost"},  // OPATMOST
    {expr::PLTERM,           prec::CALL,           "pl term"},  // OPPLTERM
    // OPIFSYM - not supported yet
    {expr::UNKNOWN,          prec::UNKNOWN,        "string if-then-else"},
    {expr::LOGICAL_COUNT,    prec::CALL,           "exactly"},  // OPEXACTLY
    {expr::LOGICAL_COUNT,    prec::CALL,           "!atleast"},  // OPNOTATLEAST
    {expr::LOGICAL_COUNT,    prec::CALL,           "!atmost"},  // OPNOTATMOST
    {expr::LOGICAL_COUNT,    prec::CALL,           "!exactly"},  // OPNOTEXACTLY
    {expr::ITERATED_LOGICAL, prec::CALL,           "forall"},  // ANDLIST
    {expr::ITERATED_LOGICAL, prec::CALL,           "exists"},  // ORLIST
    {expr::IMPLICATION,      prec::IMPLICATION,    "==>"},  // OPIMPELSE
    {expr::BINARY_LOGICAL,   prec::IFF,            "<==>"},  // OP_IFF
    {expr::ALLDIFF,          prec::CALL,           "alldiff"},  // OPALLDIFF
    {expr::BINARY,           prec::EXPONENTIATION, "^"},  // OP1POW
    {expr::UNARY,            prec::EXPONENTIATION, "^2"},  // OP2POW
    {expr::BINARY,           prec::EXPONENTIATION, "^"},  // OPCPOW
    // OPFUNCALL
    {expr::CALL,             prec::CALL,           "function call"},
    {expr::CONSTANT,         prec::PRIMARY,        "number"},  // OPNUM
    {expr::STRING,           prec::PRIMARY,        "string"},  // OPHOL
    {expr::VARIABLE,         prec::PRIMARY,        "variable"}  // OPVARVAL
};
}
