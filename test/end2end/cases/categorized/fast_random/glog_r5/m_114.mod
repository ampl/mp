var x{1..6} integer >= 0 <= 3;

minimize o: -2*x[1] - 4*x[2] + 6*x[3] + 6*x[4] + x[5] + 2*x[6];

s.t. c1: not ((x[3] > 3) and (x[5] >= 0));
