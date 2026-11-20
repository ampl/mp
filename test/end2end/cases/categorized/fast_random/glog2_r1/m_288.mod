var x{1..4} integer >= 0 <= 1;

maximize o: -2*x[1] - 3*x[2] - 2*x[3] - 2*x[4];

s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4]) > 0;
s.t. c2: (if x[3] != 1 then x[1] else x[4]) != 1;
s.t. c3: max(x[1], x[2]) - min(x[4], x[4]) > 0;
s.t. c4: numberof 0 in (x[1], x[2], x[3], x[4]) <= 1;
