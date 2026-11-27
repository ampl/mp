var x{1..5} integer >= -2 <= 0;

maximize o: 2*x[1] - 4*x[2] - x[3] - x[4] - 2*x[5];

s.t. c1: max(x[2], x[2]) - min(x[1], x[4]) != -2;
