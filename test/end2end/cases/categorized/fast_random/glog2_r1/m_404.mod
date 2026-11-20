var x{1..5} integer >= -2 <= 3;

maximize o: -5*x[1] - 2*x[2] - 3*x[3] + 0*x[4] + x[5];

s.t. c1: max(x[1], x[2]) - min(x[3], x[5]) == 1;
