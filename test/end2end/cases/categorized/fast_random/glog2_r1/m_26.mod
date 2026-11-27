var x{1..3} integer >= -2 <= 2;

maximize o: 0*x[1] + 0*x[2] + 0*x[3];

s.t. c1: alldiff(x[3], x[1]);
