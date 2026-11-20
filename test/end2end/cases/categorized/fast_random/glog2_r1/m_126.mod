var x{1..4} integer >= 0 <= 3;

minimize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4];

s.t. c1: numberof 1 in (x[1], x[2], x[3], x[4]) < 0;
