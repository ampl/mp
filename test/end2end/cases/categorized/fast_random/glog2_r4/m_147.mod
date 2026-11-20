var x{1..4} integer >= -2 <= 2;

maximize o: -3000*x[1] + 17*x[2] - 3000*x[3] - 999*x[4];

s.t. c1: ((if x[4] != 2 then 1 else 0) + (if x[1] >= 2 then 1 else 0)) >= 1;
s.t. c2: numberof 2 in (x[1], x[2], x[3], x[4]) != 2;
s.t. c3: abs(x[4] - x[3]) > 0;
s.t. c4: not ((x[1] != -1) and (x[1] == 0));
