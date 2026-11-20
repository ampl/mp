var x{1..6} integer >= 0 <= 3;

minimize o: 3*x[1] + 2*x[2] - x[3] + x[4] + 0*x[5] + 0*x[6];

s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) == 2;
s.t. c2: not ((x[2] == 0) and (x[2] != 1));
