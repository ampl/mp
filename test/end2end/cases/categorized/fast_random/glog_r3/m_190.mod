var x{1..5} integer >= -2 <= 0;
minimize o: (2)*x[1] + (4)*x[2] + (3)*x[3] + (0)*x[4] + (0)*x[5];
s.t. c1: ((if x[1] > -2 then 1 else 0) + (if x[5] >= -2 then 1 else 0)) >= 1;
s.t. c2: max(x[3], x[4]) - min(x[4], x[1]) == 0;
s.t. c3: ((if x[5] != 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
