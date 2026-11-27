var x{1..6} integer >= 0 <= 4;

minimize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5] + 0*x[6];

s.t. c1: (x[6] == 2) or (x[4] <= 0) or (x[5] <= 2);
