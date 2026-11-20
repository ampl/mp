var x{1..4} integer >= -2 <= 1;

maximize o: -5*x[1] + x[2] - 5*x[3] - 5*x[4];

s.t. c1: not ((x[3] > 1) and (x[4] >= 1));
s.t. c2: ((if x[3] == -1 then 1 else 0) + (if x[2] >= -1 then 1 else 0)) >= 1;
s.t. c3: numberof 1 in (x[1], x[2], x[3], x[4]) != 2;
s.t. c4: max(x[2], x[3]) - min(x[1], x[3]) != 1;
