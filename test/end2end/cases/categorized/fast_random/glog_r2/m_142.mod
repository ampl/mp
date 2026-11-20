var x{1..4} integer >= 0 <= 3;

maximize o: -4*x[1] - x[2] + 2*x[3] + x[4];

s.t. c1: (x[2] > 2) or (x[3] >= 3) or (x[2] < 1);
