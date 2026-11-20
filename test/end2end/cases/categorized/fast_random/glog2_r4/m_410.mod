var x{1..3} integer >= 0 <= 2;

maximize o: 4*x[1] - 3*x[2] - 4*x[3];

s.t. c1: max(x[3], x[1]) - min(x[3], x[3]) <= 1;
