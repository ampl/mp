var x{1..4} integer >= -2 <= 2;

maximize o: 4*x[1] - 5*x[2] - 2*x[3] - 5*x[4];

s.t. c1: max(x[3], x[1]) - min(x[3], x[3]) <= 3;
