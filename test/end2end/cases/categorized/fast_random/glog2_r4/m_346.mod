var x{1..4} integer >= 0 <= 4;

maximize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4];

s.t. c1: max(x[3], x[4]) - min(x[2], x[3]) != 3;
