var x{1..4} integer >= -2 <= 1;

maximize o: 3*x[1] - 5*x[2] - 5*x[3] - 2*x[4];

s.t. c1: not ((x[3] > 1) and (x[3] != 0));
