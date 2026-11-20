var x{1..6} integer >= 0 <= 1;

maximize o: -x[1] + 4*x[2] + 6*x[3] + 6*x[4] + 4*x[5] + 0*x[6];

s.t. c1: alldiff(x[5], x[4], x[3], x[1]);
