var x{1..5} integer >= -2 <= 0;

minimize o: -4*x[1] - 3*x[2] + 6*x[3] + 4*x[4] - x[5];

s.t. c1: not ((x[2] <= 0) and (x[2] != 0));
