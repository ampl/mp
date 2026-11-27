var x{1..6} integer >= 0 <= 4;

minimize o: 5*x[1] + 5*x[2] + 4*x[3] - 2*x[4] + x[5] + 0*x[6];

s.t. c1: (x[3] > 4) or (x[4] <= 1) or (x[1] < 0);
