var x{1..5} integer >= 0 <= 1;

maximize o: -4*x[1] + 3*x[2] + 2*x[3] + 6*x[4] + 4*x[5];

s.t. c1: ((if x[5] >= 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
s.t. c2: exists {i in 1..5} x[i] == 1;
