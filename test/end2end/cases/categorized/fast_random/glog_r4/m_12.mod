var x{1..6} integer >= 0 <= 3;

maximize o: 2*x[1] + x[2] + 3*x[3] - x[4] - 5*x[5] + 6*x[6];

s.t. c1: not ((x[5] >= 2) and (x[3] < 3));
