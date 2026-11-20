var x{1..3} integer >= 0 <= 4;

minimize o: -2*x[1] + x[2] - 5*x[3];

s.t. c1: (x[3] >= -1) or (x[3] > 5) or (x[2] >= 0);
