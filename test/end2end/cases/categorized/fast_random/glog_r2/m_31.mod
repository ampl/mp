var x{1..4} integer >= -2 <= 1;

maximize o: x[1] + 3*x[2] + 0*x[3] - 5*x[4];

s.t. c1: (x[3] <= 0) or (x[3] != 1) or (x[3] >= -1);
s.t. c2: numberof -2 in (x[1], x[2], x[3], x[4]) >= 2;
s.t. c3: abs(x[1] - x[1]) > 1;
s.t. c4: not ((x[3] != -1) and (x[1] == -1));
