var x{1..3} integer >= -2 <= 3;

minimize o: 5*x[1] + 0*x[2] + x[3];

s.t. c1: (if x[3] < 1 then x[2] else x[2]) > 3;
s.t. c2: numberof 2 in (x[1], x[2], x[3]) < 2;
