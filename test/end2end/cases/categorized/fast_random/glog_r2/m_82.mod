var x{1..6} integer >= -2 <= 2;

minimize o: 3*x[1] + 5*x[2] + 5*x[3] + 2*x[4] - 5*x[5] - x[6];

s.t. c1: not ((x[5] == 2) and (x[1] == 1));
