var x{1..6} integer >= 0 <= 3;
minimize o: (-2)*x[1] + (-3)*x[2] + (-1)*x[3] + (0)*x[4] + (5)*x[5] + (5)*x[6];
s.t. c1: ((if x[3] == 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
s.t. c2: x[2] == 3;
