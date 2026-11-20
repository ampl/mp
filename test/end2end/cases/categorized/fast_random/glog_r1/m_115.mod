var x{1..3} integer >= 0 <= 4;

minimize o: x[1] + x[2] - 4*x[3];

s.t. c1: numberof 2 in (x[1], x[2], x[3]) >= 2;
