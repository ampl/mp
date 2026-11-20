var x{1..3} integer >= 0 <= 2;

minimize o: -4*x[1] - 5*x[2] - x[3];

s.t. c1: abs(x[2] - x[1]) == 0;
s.t. c2: max(x[2], x[2]) - min(x[1], x[3]) < 2;
s.t. c3: numberof 1 in (x[1], x[2], x[3]) >= 2;
s.t. c4: max(x[1], x[3]) - min(x[3], x[3]) < 1;
