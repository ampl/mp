# int_ne_06_ctx_redef.mod
#
# Test that x<=5 is redefined after all
# its contexts are known

var x >=-3 <=11 integer;
var b binary;
var c binary;

minimize Obj: x + b + c;

s.t. ConIMPL_01: !(x<=5) or c;

s.t. ConNE_01: x != 6 <==> b;