var x{1..5} integer >= -3 <= 1;

minimize o: -5*x[1] - 2*x[2] + x[3] - 5*x[4] - 5*x[5];

s.t. c1: exists {i in 1..5} x[i] == -2;
s.t. c2: exists {i in 1..5} x[i] == 1;
s.t. c3: ((if x[4] < 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
