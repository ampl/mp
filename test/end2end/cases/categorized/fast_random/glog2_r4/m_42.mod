var x{1..6} integer >= 0 <= 3;

minimize o: x[1] + 6*x[2] - 4*x[3] - 3*x[4] + 4*x[5] - x[6];

s.t. c1: max(x[2], x[4]) - min(x[5], x[1]) < 3;
