var x{1..3} integer >= 0 <= 3;

minimize o: 5*x[1] + 4*x[2] - 4*x[3];

s.t. c1: alldiff(x[1], x[2], x[3]);
s.t. c2: not ((x[1] <= 3) and (x[1] >= 2));
s.t. c3: alldiff(x[2], x[3], x[1]);
