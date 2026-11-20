var x{1..5} integer >= -2 <= 1;

maximize o: 6*x[1] - 5*x[2] + 5*x[3] + 6*x[4] + 4*x[5];

s.t. c1: max(x[4], x[4]) - min(x[4], x[2]) > 0;
