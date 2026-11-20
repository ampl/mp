var x{1..3} integer >= -2 <= 1;

maximize o: 5*x[1] + 5*x[2] + 5*x[3];

s.t. c1: max(x[2], x[3]) - min(x[2], x[1]) > 1;
