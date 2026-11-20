var x{1..7} integer >= 0 <= 1;

minimize o: 0*x[1] - 4*x[2] - 4*x[3] - 4*x[4] - 2*x[5] + 2*x[6] + 2*x[7];

s.t. c1: forall {i in 1..6} x[i] != x[i+1];
s.t. c2: numberof 1 in (x[1], x[2], x[3], x[4], x[5], x[6], x[7]) == 3;
