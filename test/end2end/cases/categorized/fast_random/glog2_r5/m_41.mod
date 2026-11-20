var x{1..6} integer >= -2 <= 0;

maximize o: -x[1] + 5*x[2] + 4*x[3] - 3*x[4] - x[5] + 5*x[6];

s.t. c1: alldiff(x[2], x[6], x[5]);
s.t. c2: alldiff(x[5], x[2], x[3]);
s.t. c3: not ((x[4] == 0) and (x[5] == -2));
