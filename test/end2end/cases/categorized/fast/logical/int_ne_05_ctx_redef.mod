# int_ne_05_ctx_redef.mod
#
# Test that x<=5 is redefined after all
# its contexts are known #248

var x >=-3 <=11 integer;
var b binary;
var b1: binary;

minimize Obj: -10*b - 10*b1 + if x<=5 then 1;

s.t. ConNE_01: x != 6 <==> b;

s.t. ConImpl_01: b1 ==> x <= 6;
