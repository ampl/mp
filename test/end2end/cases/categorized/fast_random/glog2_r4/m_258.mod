var x{1..4} integer >= 0 <= 1;

maximize o: -2*x[1] - 2*x[2] - 2*x[3] + 3*x[4];

s.t. c1: not ((x[3] > 0) and (x[3] == 1));
