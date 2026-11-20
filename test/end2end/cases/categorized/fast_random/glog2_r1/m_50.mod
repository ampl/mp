var x{1..6} integer >= 0 <= 3;

maximize o: -x[1] + 3*x[2] - 5*x[3] - 5*x[4] + 0*x[5] + 0*x[6];

s.t. c1: not ((x[1] == 2) and (x[4] == 3));
