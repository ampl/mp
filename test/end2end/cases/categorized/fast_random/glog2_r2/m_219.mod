var x{1..3} integer >= 0 <= 4;

minimize o: -x[1] - 5*x[2] - 3*x[3];

s.t. c1: alldiff(x[1], x[3]);
s.t. c2: alldiff(x[1], x[3]);
