var x{1..3} integer >= -3 <= 2;

maximize o: -5*x[1] - 5*x[2] - x[3];

s.t. c1: alldiff(x[2], x[1]);
