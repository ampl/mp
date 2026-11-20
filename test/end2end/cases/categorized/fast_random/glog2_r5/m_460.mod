var x{1..6} integer >= 0 <= 4;

minimize o: 0*x[1] + 2*x[2] + 3*x[3] + 6*x[4] + 4*x[5] + 0*x[6];

s.t. c1: max(x[6], x[2]) - min(x[4], x[5]) != 4;
s.t. c2: abs(x[6] - x[2]) == 0;
