var x{1..5} integer >= 0 <= 3;
minimize o: (-2)*x[1] + (6)*x[2] + (2)*x[3] + (1)*x[4] + (1)*x[5];
s.t. c1: ((if x[3] < 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
