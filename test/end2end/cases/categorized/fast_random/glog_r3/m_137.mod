var x{1..3} integer >= 0 <= 4;

maximize o: -3*x[1] - 5*x[2] - 5*x[3];

s.t. c1: alldiff(x[2], x[3], x[1]);
