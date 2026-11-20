var x{1..6} integer >= -2 <= 0;

minimize o: 2*x[1] - 2*x[2] - 4*x[3] + 4*x[4] + 3*x[5] - 3*x[6];

s.t. c1: alldiff(x[6], x[1], x[2], x[3]);
