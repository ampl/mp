# -------------------------------------------------------------
# funcbnd_log_02_int.mod: log of an integer variable, approximated
#   piecewise-linearly. With exact result bounds (cvt:pre:funcbndtol=0)
#   the domain clipped by the function values ends at
#   exp(log(5)) = 4.999...: the breakpoint x = 5 must not be lost.
#   Optimum x = 4, objective log(4) - 1.08 = 0.3062943611.
# -------------------------------------------------------------

var x integer >= 1, <= 5;

maximize Obj: log(x) - 0.27*x;
