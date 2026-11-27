var x{1..4} integer >= 0 <= 1;

maximize o: -5*x[1] + x[2] + 2*x[3] - 2*x[4];

s.t. c1: max(x[1], x[2]) - min(x[3], x[3]) >= 1;
