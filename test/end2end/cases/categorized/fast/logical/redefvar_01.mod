######################################################
## redefvar_01.mod:
## test redefinitions of auxiliary variables.
## As of v4.0.4, it seems to work,
## because MakeComplementVar() uses Not()
## andLinearFuncCon's are inlined.
######################################################

var x >=0 <=10;
var y >=0 <=10;
var b binary;

s.t. DefB: b <==> x >= 5;

s.t. NotBRedefined: y <= 4 || if 1-b then 1;

s.t. NotB: y >= 6 || not(b);
