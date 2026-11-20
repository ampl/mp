var x{1..3} integer >= 0 <= 2;

maximize o: -4*x[1] + x[2] - 4*x[3];

s.t. c1: not ((x[2] < 0) and (x[3] != 0));
s.t. c2: not ((x[2] < 0) and (x[3] != 0));
