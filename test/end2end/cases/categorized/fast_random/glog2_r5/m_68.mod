var x{1..5} integer >= 0 <= 1;

maximize o: 2*x[1] - 4*x[2] - 4*x[3] + 4*x[4] - 5*x[5];

s.t. c1: exists {i in 1..5} x[i] == 1;
s.t. c2: (x[5] != 1) or (x[5] > 0) or (x[1] > 0);
s.t. c3: ((if x[2] != 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
s.t. c4: numberof 1 in (x[1], x[2], x[3], x[4], x[5]) >= 2;
