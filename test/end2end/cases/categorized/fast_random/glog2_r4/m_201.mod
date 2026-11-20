var x{1..4} integer >= 0 <= 2;

maximize o: -5*x[1] - 5*x[2] - 3*x[3] - 2*x[4];

s.t. c1: alldiff(x[1], x[4], x[2], x[3]);
