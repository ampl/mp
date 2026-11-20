var x{1..6} integer >= -2 <= 1;
minimize o: (5)*x[1] + (0)*x[2] + (-3)*x[3] + (4)*x[4] + (-3)*x[5] + (3)*x[6];
s.t. c1: abs(x[5] - x[4]) > 2;
s.t. c2: ((if x[5] == -2 then 1 else 0) + (if x[2] >= -2 then 1 else 0)) >= 1;
s.t. c3: (if x[6] < 1 then x[1] else x[3]) < 1;
