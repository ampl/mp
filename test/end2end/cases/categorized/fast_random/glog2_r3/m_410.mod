var x{1..3} integer >= 0 <= 3;

minimize o: -3*x[1] + 2*x[2] - 3*x[3];

s.t. c1: max(x[2], x[1]) - min(x[1], x[3]) != 0;
s.t. c2: numberof 3 in (x[1], x[2], x[3]) >= 2;
