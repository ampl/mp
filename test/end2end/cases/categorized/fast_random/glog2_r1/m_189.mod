var x{1..3} integer >= -2 <= 2;

minimize o: -x[1] + 4*x[2] - 4*x[3];

s.t. c1: numberof -1 in (x[1], x[2], x[3]) >= 0;
s.t. c2: numberof -1 in (x[1], x[2], x[3]) >= 0;
