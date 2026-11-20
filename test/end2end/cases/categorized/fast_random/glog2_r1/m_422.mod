var x{1..6} integer >= -2 <= 2;

minimize o: 2*x[1] + 3*x[2] + x[3] + 3*x[4] - 5*x[5] + 3*x[6];

s.t. c1: (x[5] < -1) or (x[4] > -2) or (x[5] >= 2);
