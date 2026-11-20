var x{1..4} integer >= 0 <= 2;

maximize o: 0*x[1] + 3*x[2] + 4*x[3] + 6*x[4];

s.t. c1: max(x[4], x[4]) - min(x[4], x[3]) >= 2;
