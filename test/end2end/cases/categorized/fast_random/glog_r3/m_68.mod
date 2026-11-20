var x{1..5} integer >= -2 <= 0;

minimize o: x[1] + 2*x[2] + x[3] - 2*x[4] - 3*x[5];

s.t. c1: max(x[4], x[3]) - min(x[2], x[3]) != 0;
