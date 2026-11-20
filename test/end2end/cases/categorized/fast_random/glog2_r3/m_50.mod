var x{1..6} integer >= 0 <= 3;

minimize o: 5*x[1] + 0*x[2] - 2*x[3] + 0*x[4] + 4*x[5] + 3*x[6];

s.t. c1: (x[2] >= 1) or (x[1] > 0) or (x[6] <= 3);
