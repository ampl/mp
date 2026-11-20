var x{1..5} integer >= -2 <= 2;

maximize o: -4*x[1] - 2*x[2] + x[3] + x[4] + 0*x[5];

s.t. c1: not ((x[5] != -2) and (x[5] != 1));
