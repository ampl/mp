var x{1..6} integer >= -2 <= 2;

minimize o: 3*x[1] - 4*x[2] + 6*x[3] + 0*x[4] + x[5] - 5*x[6];

s.t. c1: alldiff(x[2], x[3], x[1], x[5]);
