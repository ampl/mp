var x{1..6} integer >= -2 <= 1;

minimize o: 6*x[1] + 4*x[2] + 0*x[3] - 5*x[4] + 0*x[5] + 4*x[6];

s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) == 2;
