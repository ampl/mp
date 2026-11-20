var x{1..5} integer >= -2 <= 0;

minimize o: -3*x[1] + x[2] + x[3] - 3*x[4] - 3*x[5];

s.t. c1: numberof -1 in (x[1], x[2], x[3], x[4], x[5]) < 0;
