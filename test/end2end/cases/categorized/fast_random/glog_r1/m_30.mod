var x{1..6} integer >= 0 <= 2;

minimize o: 6*x[1] + 2*x[2] + 3*x[3] - 4*x[4] - 4*x[5] + 4*x[6];

s.t. c1: alldiff(x[5], x[4]);
