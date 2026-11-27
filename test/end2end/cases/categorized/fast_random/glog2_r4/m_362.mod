var x{1..3} integer >= -2 <= 0;

maximize o: -3*x[1] - 5*x[2] + 4*x[3];

s.t. c1: not ((x[2] > -1) and (x[2] >= 0));
s.t. c2: alldiff(x[1], x[3]);
