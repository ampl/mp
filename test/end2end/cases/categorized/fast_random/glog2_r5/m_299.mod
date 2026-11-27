var x{1..3} integer >= 0 <= 2;

maximize o: x[1] + 2*x[2] + 3*x[3];

s.t. c1: (x[1] != 1) or (x[3] >= 0) or (x[1] >= 0);
s.t. c2: (x[1] != 1) or (x[3] >= 0) or (x[1] >= 0);
