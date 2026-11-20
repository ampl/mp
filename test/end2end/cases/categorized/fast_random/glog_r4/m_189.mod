var x{1..3} integer >= -2 <= 0;

maximize o: -3*x[1] - 4*x[2] + 5*x[3];

s.t. c1: max(x[2], x[2]) - min(x[3], x[3]) < 0;
