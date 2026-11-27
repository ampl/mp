var x{1..5} integer >= -1 <= 1;

maximize o: 3*x[1] + 5*x[2] + 3*x[3] - 4*x[4] - 5*x[5];

s.t. c1: not (x[4] != 0);
s.t. c2: (x[4] != 0) or (x[3] <= -1);
