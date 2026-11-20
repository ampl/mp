var x{1..3} integer >= 0 <= 2;

maximize o: -4*x[1] + 3*x[2] + 3*x[3];

s.t. c1: not ((x[2] == 2) and (x[3] <= 1));
s.t. c2: (x[1] > 1) or (x[3] != 2) or (x[3] >= 1);
