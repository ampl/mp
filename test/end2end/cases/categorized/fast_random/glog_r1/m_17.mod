var x{1..4} integer >= 0 <= 4;

maximize o: 6*x[1] + 0*x[2] + x[3] - 3*x[4];

s.t. c1: (x[2] != 4) or (x[1] <= 0) or (x[2] > 0);
