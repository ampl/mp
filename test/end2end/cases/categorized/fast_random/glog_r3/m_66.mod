var x{1..3} integer >= -2 <= 1;

minimize o: 4*x[1] + 5*x[2] + x[3];

s.t. c1: numberof 0 in (x[1], x[2], x[3]) >= 1;
s.t. c2: not ((x[2] == 0) and (x[1] != -1));
