var x{1..6} integer >= 0 <= 1;

minimize o: -2*x[1] + 4*x[2] + 2*x[3] + 0*x[4] + 5*x[5] - x[6];

s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) > 1;
s.t. c2: (x[2] == 0) or (x[5] > 1) or (x[2] > 0);
