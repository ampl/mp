var x{1..4} integer >= 0 <= 4;

maximize o: 2*x[1] - 5*x[2] - 2*x[3] + 0*x[4];

s.t. c1: max(x[2], x[2]) - min(x[4], x[4]) <= 0;
s.t. c2: (x[4] > 2) or (x[1] >= 0) or (x[3] != 0);
s.t. c3: numberof 4 in (x[1], x[2], x[3], x[4]) != 3;
s.t. c4: ((if x[4] != 3 then 1 else 0) + (if x[3] >= 3 then 1 else 0)) >= 1;
