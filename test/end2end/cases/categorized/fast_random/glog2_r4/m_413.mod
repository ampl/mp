var x{1..4} integer >= -2 <= 0;

minimize o: x[1] + 0*x[2] - x[3] - 2*x[4];

s.t. c1: numberof -1 in (x[1], x[2], x[3], x[4]) <= 1;
s.t. c2: max(x[4], x[3]) - min(x[3], x[2]) >= 0;
s.t. c3: x[4] == -1;
s.t. c4: x[1] == 0;
