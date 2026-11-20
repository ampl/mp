var x{1..6} integer >= 0 <= 3;

minimize o: x[1] - 4*x[2] + 6*x[3] + x[4] + 0*x[5] + 5*x[6];

s.t. c1: max(x[1], x[4]) - min(x[5], x[4]) < 1;
