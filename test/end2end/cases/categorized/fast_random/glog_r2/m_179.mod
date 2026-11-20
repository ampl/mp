var x{1..6} integer >= 0 <= 3;

minimize o: -3*x[1] + x[2] + 0*x[3] + x[4] - 3*x[5] + 3*x[6];

s.t. c1: not ((x[6] != 1) and (x[1] >= 2));
