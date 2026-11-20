var x{1..6} integer >= 0 <= 2;
minimize o: (1)*x[1] + (2)*x[2] + (5)*x[3] + (0)*x[4] + (-4)*x[5] + (-3)*x[6];
s.t. c1: not ((x[6] < 1) and (x[2] > 2));
s.t. c2: abs(x[4] - x[5]) <= 2;
s.t. c3: ((if x[1] < 2 then 1 else 0) + (if x[1] >= 2 then 1 else 0)) >= 1;
