var x{1..4} integer >= 0 <= 3;

maximize o: -4*x[1] + 0*x[2] + 5*x[3] + 3*x[4];

s.t. c1: max(x[1], x[4]) - min(x[4], x[3]) <= 1;
s.t. c2: numberof 2 in (x[1], x[2], x[3], x[4]) != 2;
s.t. c3: ((if x[2] > 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c4: exists {i in 1..4} x[i] == 0;
