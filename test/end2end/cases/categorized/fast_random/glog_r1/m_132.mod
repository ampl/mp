var x{1..4} integer >= 0 <= 4;

minimize o: -3*x[1] + 5*x[2] - 4*x[3] + 2*x[4];

s.t. c1: not ((x[2] != 1) and (x[4] <= 2));
