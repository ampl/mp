var x{1..5} integer >= -2 <= 0;

maximize o: 6*x[1] - x[2] + 4*x[3] + 5*x[4] + x[5];

s.t. c1: not ((x[5] == 0) and (x[3] >= 0));
s.t. c2: not ((x[5] == 0) and (x[3] >= 0));
