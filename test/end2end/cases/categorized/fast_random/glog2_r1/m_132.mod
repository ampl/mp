var x{1..3} integer >= 0 <= 3;

minimize o: -2*x[1] + 6*x[2] + x[3];

s.t. c1: numberof 1 in (x[1], x[2], x[3]) < 0;
s.t. c2: abs(x[2] - x[1]) < 0;
