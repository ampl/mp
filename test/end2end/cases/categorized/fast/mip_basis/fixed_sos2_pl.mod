# -------------------------------------------------------------
# fixed_sos2_pl.mod
# Nonconvex piecewise-linear terms that AMPL linearizes with SOS2
# (option pl_linearize 1), plus logical constraints. Optimum 0.
# The fixed LP solved for duals (mip:basis, cplexasl mipbasis)
# has no SOS constraints; its optimum -3.3 at x = 2 is that of the
# convex relaxation of the terms, where obj evaluates to 2.7.
# -------------------------------------------------------------
var x >= -3 <= 5;
var y >= -2 <= 4;
var b binary;

minimize obj:
   <<-1, 1, 2; -1, 0.5, 2, -1>> x + 3*b + 0.1*x;

s.t. c1: 2*y + <<0; 1, -2>> y >= -1;
s.t. c2: <<-1, 1, 2; -1, 0.5, 2, -1>> x + y <= 4;

s.t. lc: b = 1 ==> x >= 2 else y <= 1;
s.t. lc2: b = 0 ==> x + y >= 1;
