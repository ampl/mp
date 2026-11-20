var x{1..6} integer >= 0 <= 3;

minimize o: 3*x[1] - 2*x[2] - x[3] - 5*x[4] + 2*x[5] - 5*x[6];

s.t. c1: max(x[6], x[6]) - min(x[6], x[1]) != 3;
