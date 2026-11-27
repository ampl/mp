var x{1..5} integer >= 0 <= 3;

maximize o: 6*x[1] + 2*x[2] + 0*x[3] + 2*x[4] + 2*x[5];

s.t. c1: alldiff(x[1], x[2], x[4], x[5]);
