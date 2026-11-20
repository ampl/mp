var x{1..6} integer >= 0 <= 1;

minimize o: -5*x[1] + 5*x[2] - 5*x[3] - 5*x[4] + 3*x[5] + 6*x[6];

s.t. c1: alldiff(x[4], x[5], x[6]);
