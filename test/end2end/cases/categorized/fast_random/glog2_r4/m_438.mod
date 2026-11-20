var x{1..7} integer >= 0 <= 1;

minimize o: -3*x[1] + 0*x[2] + 0*x[3] + x[4] - x[5] + 2*x[6] - 5*x[7];

s.t. c1: abs(x[1] - x[5]) <= 1;
s.t. c2: max(x[2], x[1]) - min(x[5], x[2]) < 1;
s.t. c3: forall {i in 1..6} x[i] != x[i+1];
