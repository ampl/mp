var x{1..6} integer >= 0 <= 4;

maximize o: 0*x[1] - 3*x[2] - 2*x[3] + 2*x[4] - 4*x[5] - 5*x[6];

s.t. c1: not ((x[6] > 0) and (x[1] != 3));
