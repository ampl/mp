var x{1..6} integer >= 0 <= 4;

maximize o: 4*x[1] - 2*x[2] - x[3] - 3*x[4] - x[5] + x[6];

s.t. c1: not ((x[2] != 3) and (x[5] != 0));
