var x{1..3} integer >= 0 <= 3;

minimize o: 4*x[1] - 5*x[2] + 0*x[3];

s.t. c1: abs(x[1] - x[1]) != 1;
s.t. c2: numberof 1 in (x[1], x[2], x[3]) <= 3;
s.t. c3: numberof 0 in (x[1], x[2], x[3]) >= 1;
s.t. c4: abs(x[2] - x[1]) > 1;
