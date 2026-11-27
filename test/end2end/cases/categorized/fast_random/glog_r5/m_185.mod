var x{1..3} integer >= 0 <= 2;

maximize o: -3*x[1] - 3*x[2] + 4*x[3];

s.t. c1: alldiff(x[2], x[1], x[3]);
