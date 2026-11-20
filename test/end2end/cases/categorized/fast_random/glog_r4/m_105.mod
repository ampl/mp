var x{1..6} integer >= 0 <= 4;
minimize o: (1)*x[1] + (2)*x[2] + (5)*x[3] + (6)*x[4] + (0)*x[5] + (1)*x[6];
s.t. c1: ((if x[2] != 2 then 1 else 0) + (if x[1] >= 2 then 1 else 0)) >= 1;
s.t. c2: max(x[2], x[3]) - min(x[3], x[6]) == 4;
s.t. c3: abs(x[4] - x[5]) >= 1;
s.t. c4: max(x[5], x[5]) - min(x[1], x[6]) <= 1;
