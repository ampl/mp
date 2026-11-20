var x{1..5} integer >= 0 <= 3;

maximize o: -2*x[1] - 5*x[2] + 0*x[3] - 2*x[4] + 5*x[5];

s.t. c1: alldiff(x[3], x[1], x[2], x[5]);
