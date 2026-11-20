var x{1..3} integer >= 0 <= 4;

minimize o: -4*x[1] + x[2] + 0*x[3];

s.t. c1: numberof 2 in (x[1], x[2], x[3]) > 0;
s.t. c2: exists {i in 1..3} x[i] == 1;
