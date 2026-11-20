var x{1..5} integer >= -3 <= 2;

minimize o: -x[1] - 5*x[2] + 2*x[3] - 5*x[4] + 6*x[5];

s.t. c1: max(x[3], x[4]) - min(x[3], x[3]) == 1;
