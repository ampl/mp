var x{1..5} integer >= -2 <= 0;

minimize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5];

s.t. c1: ((if x[5] <= -1 then 1 else 0) + (if x[2] >= -1 then 1 else 0)) >= 1;
s.t. c2: exists {i in 1..5} x[i] == -1;
s.t. c3: max(x[2], x[5]) - min(x[5], x[4]) <= -1;
