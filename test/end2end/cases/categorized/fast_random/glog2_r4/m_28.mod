var x{1..4} integer >= 0 <= 1;

minimize o: -2*x[1] + 6*x[2] + 0*x[3] - 4*x[4];

s.t. c1: not ((x[4] == 0) and (x[3] == 1));
