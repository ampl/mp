var x{1..6} integer >= -2 <= 1;

minimize o: x[1] + 6*x[2] - 5*x[3] - 5*x[4] + 6*x[5] + 0*x[6];

s.t. c1: (x[4] >= 2) or (x[1] >= -2) or (x[4] < -2);
