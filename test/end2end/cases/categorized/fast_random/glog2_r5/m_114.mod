var x{1..5} integer >= -2 <= 3;

maximize o: x[1] + x[2] + 0*x[3] - x[4] - 5*x[5];

s.t. c1: max(x[4], x[5]) - min(x[5], x[3]) < 2;
