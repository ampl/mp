var x{1..5} integer >= 0 <= 2;

minimize o: 6*x[1] + 5*x[2] + 2*x[3] + x[4] + 3*x[5];

s.t. c1: alldiff(x[4], x[2], x[1]);
s.t. c2: alldiff(x[4], x[2], x[1]);
