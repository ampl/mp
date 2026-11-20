var x{1..4} integer >= -3 <= 2;

minimize o: -x[1] + x[2] + 6*x[3] - 2*x[4];

s.t. c1: alldiff(x[3], x[1]);
