var x{1..4} integer >= 0 <= 3;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4];
s.t. c1: ((if x[4] > 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
s.t. c2: max(x[4], x[4]) - min(x[1], x[4]) >= 1;
s.t. c3: (x[3] >= 0) or (x[1] > 3) or (x[4] > 3);
