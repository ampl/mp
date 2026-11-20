var x{1..6} integer >= -2 <= 2;

minimize o: -2*x[1] + 5*x[2] + x[3] + 5*x[4] + x[5] - x[6];

s.t. c1: alldiff(x[5], x[6]);
s.t. c2: alldiff(x[5], x[6]);
