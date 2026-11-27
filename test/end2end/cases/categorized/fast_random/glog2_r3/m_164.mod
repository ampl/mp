var x{1..3} integer >= -2 <= 3;

minimize o: 0*x[1] + 0*x[2] - 3*x[3];

s.t. c1: (x[1] >= 0) or (x[2] <= -2) or (x[3] != 3);
