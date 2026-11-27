var x{1..6} integer >= 0 <= 4;

maximize o: 3*x[1] - 4*x[2] - 4*x[3] + 4*x[4] + x[5] + 0*x[6];

s.t. c1: alldiff(x[5], x[6], x[1], x[2]);
