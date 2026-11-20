var x{1..6} integer >= 0 <= 1;

maximize o: x[1] - 5*x[2] - 4*x[3] - x[4] - x[5] + 3*x[6];

s.t. c1: alldiff(x[2], x[3], x[1], x[4]);
s.t. c2: alldiff(x[3], x[5], x[6], x[4]);
