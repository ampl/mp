var x{1..4} integer >= -2 <= 3;

maximize o: -4*x[1] + 3*x[2] + 0*x[3] - 5*x[4];

s.t. c1: ((if x[1] <= 2 then 1 else 0) + (if x[4] >= 2 then 1 else 0)) >= 1;
s.t. c2: ((if x[2] <= 1 then 1 else 0) + (if x[4] >= 1 then 1 else 0)) >= 1;
s.t. c3: exists {i in 1..4} x[i] == 2;
s.t. c4: numberof 3 in (x[1], x[2], x[3], x[4]) > 1;
