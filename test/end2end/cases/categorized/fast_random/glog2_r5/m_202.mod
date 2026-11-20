var x{1..3} integer >= 0 <= 4;

minimize o: 3*x[1] + 3*x[2] + 4*x[3];

s.t. c1: numberof 0 in (x[1], x[2], x[3]) > 0;
s.t. c2: max(x[1], x[2]) - min(x[2], x[2]) <= 4;
s.t. c3: not ((x[3] <= 1) and (x[3] <= 0));
s.t. c4: abs(x[2] - x[2]) > 4;
