######################################################
## redefvar_02.mod:
## test redefinitions of auxiliary variables.
## As of v4.0.4, it fails, because b1+b2 is known
## before count(b1, b2) which is redefined into it.
## Needs acc:linfunccon=0 ???
######################################################

var x >=0 <=10;
var y >=0 <=10;
var b1 binary;
var b2 binary;

s.t. DefB1: b1 <==> x >= 5;
s.t. DefB2: b2 <==> x >= 7;

## For Gurobi 12, max() requires outlining b1+b2.
## This is also enforced by acc:linfunccon=0.
## This converts "aux_var = b1+b2" into the corr.
## linear constraint and it's important to have
## the correct aux_var
s.t. UseLFC: y >= 6 || max(b1+b2, x-y) >= 2.3;

s.t. UseCount: y <= 4 || count(b1, b2) <= 1;

maximize OO: y+b1+b2;

