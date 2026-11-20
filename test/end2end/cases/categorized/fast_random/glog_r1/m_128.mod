var x{1..5} integer >= 0 <= 2;

minimize o: -5*x[1] + 6*x[2] - 3*x[3] - 4*x[4] + 4*x[5];

s.t. c1: exists {i in 1..5} x[i] == 2;
s.t. c2: ((if x[3] >= 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
