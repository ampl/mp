var x{1..5} integer >= 0 <= 2;

minimize o: 0*x[1] + 2*x[2] + x[3] + 6*x[4] - 3*x[5];

s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) <= 3;
s.t. c2: x[1] == 0;
