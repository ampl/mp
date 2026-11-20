var x{1..5} integer >= -2 <= 0;

minimize o: -2*x[1] - 5*x[2] - x[3] - 3*x[4] - 3*x[5];

s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) < 0;
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) < 0;
