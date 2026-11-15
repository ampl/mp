# Duals of the fixed model of a nonlinear MIP (mip:basis=1): the tight row
# c has dual 1, the slack row l has 0. With the NLP presolve left on in the
# fixed solve, a driver returned them in the presolved row order (swapped).
var x >= -5, <= 5;
var z >= -5, <= 5;
var b binary;
var w integer >= 0, <= 3;
minimize o: exp(x) + exp(z) + w + 2*b;
s.t. c: x + z + w >= 2.5;
s.t. l: x + z <= 4;
s.t. xb: x >= 0.5;
