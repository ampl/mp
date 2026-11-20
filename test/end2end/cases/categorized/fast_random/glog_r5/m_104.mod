var x{1..5} integer >= 0 <= 3;
minimize o: (6)*x[1] + (-3)*x[2] + (0)*x[3] + (4)*x[4] + (5)*x[5];
s.t. c1: (x[2] < 3) or (x[3] >= 1) or (x[1] < 1);
s.t. c2: ((if x[5] > 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
