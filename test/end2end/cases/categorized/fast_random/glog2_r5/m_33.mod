var x{1..6} integer >= 0 <= 2;
minimize o: (2)*x[1] + (6)*x[2] + (-3)*x[3] + (-2)*x[4] + (1)*x[5] + (5)*x[6];
s.t. c1: ((if x[1] > 0 then 1 else 0) + (if x[6] >= 0 then 1 else 0)) >= 1;
s.t. c2: x[5] == 2;
