var x{1..4} integer >= -3 <= 1;

minimize o: 4*x[1] - 3*x[2] + 0*x[3] - 5*x[4];

s.t. c1: not ((x[4] >= -2) and (x[2] > -1));
