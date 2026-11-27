var x{1..5} integer >= -2 <= 0;

minimize o: -4*x[1] + 0*x[2] + 5*x[3] + 0*x[4] - 5*x[5];

s.t. c1: (x[5] < -1) or (x[4] != -1) or (x[5] < -2);
