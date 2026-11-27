# -------------------------------------------------------------
# cmp_tol_strict_eps_01.mod
# Strict comparison of an expression whose range (6e-6) is
# smaller than the gap cvt:mip:eps (1e-4): 1e-6*x < 3e-6 must hold,
# but the encoding requires 1e-6*x <= 3e-6 - 1e-4, which is
# infeasible. A warning names cvt:mip:eps (strict_cmp_eps).
# With cvt:mip:eps=1e-7: optimum 2.9.
# -------------------------------------------------------------
var x >= 0 <= 6;
var b binary;

maximize o: x;

s.t. d: b <==> (1e-6*x < 3e-6);
s.t. f: b >= 1;
