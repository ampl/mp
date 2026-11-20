var x{1..5} integer >= 0 <= 1;

maximize o: x[1] + 4*x[2] - 2*x[3] + 6*x[4] - 4*x[5];

s.t. c1: (x[5] < 0) or (x[2] < 1) or (x[3] >= 0);
