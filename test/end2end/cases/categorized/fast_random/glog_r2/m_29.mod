var x{1..5} integer >= 0 <= 3;

minimize o: -x[1] + 0*x[2] + 5*x[3] + 3*x[4] + x[5];

s.t. c1: abs(x[3] - x[1]) == 0;
s.t. c2: max(x[4], x[5]) - min(x[4], x[4]) > 3;
