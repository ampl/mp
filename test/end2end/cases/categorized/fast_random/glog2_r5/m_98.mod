var x{1..7} integer >= 0 <= 1;

minimize o: -5*x[1] + 4*x[2] + x[3] + 5*x[4] + 6*x[5] + x[6] + 4*x[7];

s.t. c1: not ((x[6] < 1) and (x[1] <= 0));
