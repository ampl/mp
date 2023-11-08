#ifndef NLSOL_EX_C_MODEL_H
#define NLSOL_EX_C_MODEL_H

// #include "api/c/nl-header.h"
//#include "mp/nl-opcodes.h"

/**
 * A linear model for NL Writer C API example.
 * Illustrates NL variable order (continuous -> integer),
 * linear constraints.
 *

## In AMPL, set option presolve 0; to reproduce
##
## To write NL file and name files in AMPL, use commands
##    ampl: option nl_comments 1;
##    ampl: option auxfiles rc;
##    ampl: write gmodel;

var x >=0, integer;
var y >=-17, <=504;
maximize TotalSum:
    x + 13*y;
subj to C2:
       3700*x + 0.6*y <= 3e4;
subj to C3:
       22*x + 14536*y <= 3e5;

## Initial guess
let x := 1.5;
let y := 0.11;

## A suffix
suffix zork;
let C2.zork := 5;

 *
 */


#endif // NLSOL_EX_C_MODEL_H
