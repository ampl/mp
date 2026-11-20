var x{1..4} integer >= 0 <= 3;

minimize o: 0*x[1] + 3*x[2] + 5*x[3] - 4*x[4];

s.t. c1: (x[1] < 0) or (x[3] >= 1) or (x[3] <= 2);
