var x{1..6} integer >= 0 <= 2;

minimize o: 0*x[1] + 3*x[2] + 2*x[3] - 3*x[4] + 4*x[5] + 4*x[6];

s.t. c1: not ((x[3] == 0) and (x[5] == 2));
