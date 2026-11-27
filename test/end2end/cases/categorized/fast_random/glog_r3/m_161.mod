var x{1..6} integer >= 0 <= 2;

minimize o: x[1] + 4*x[2] - x[3] - 2*x[4] + x[5] - 4*x[6];

s.t. c1: ((if x[5] <= 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
s.t. c2: exists {i in 1..6} x[i] == 0;
