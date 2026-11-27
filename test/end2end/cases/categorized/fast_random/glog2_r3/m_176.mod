var x{1..6} integer >= -2 <= 0;

minimize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5] + 0*x[6];

s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) <= 1;
s.t. c2: max(x[5], x[1]) - min(x[1], x[3]) < -2;
