var x{1..3} integer >= 0 <= 3;

maximize o: x[1] + 5*x[2] + 5*x[3];

s.t. c1: max(x[3], x[2]) - min(x[2], x[1]) != 0;
