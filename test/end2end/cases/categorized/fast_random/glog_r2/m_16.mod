var x{1..4} integer >= 0 <= 2;

minimize o: 4*x[1] + 6*x[2] + 3*x[3] + 4*x[4];

s.t. c1: alldiff(x[2], x[1], x[4], x[3]);
