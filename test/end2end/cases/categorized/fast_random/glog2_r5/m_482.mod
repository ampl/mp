var x{1..5} integer >= 0 <= 3;

maximize o: -2*x[1] + 3*x[2] + 0*x[3] - x[4] + 0*x[5];

s.t. c1: (x[4] < 2) or (x[1] != 1) or (x[1] >= 0);
