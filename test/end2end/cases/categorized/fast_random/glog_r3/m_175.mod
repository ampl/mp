var x{1..6} integer >= -2 <= 2;

minimize o: -5*x[1] + 0*x[2] + 2*x[3] + 4*x[4] + 4*x[5] - 3*x[6];

s.t. c1: not ((x[4] != -2) and (x[1] == -1));
