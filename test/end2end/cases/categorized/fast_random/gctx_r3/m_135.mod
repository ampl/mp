var x{1..4} integer >= -1 <= 1;
minimize o: (5)*x[1] + (1)*x[2] + (0)*x[3] + (-3)*x[4] + (4)*(if x[1] != x[2] then 1 else 0);
s.t. c1: not ((x[3] != 1) and (x[3] != 1));
s.t. c2: (if x[1] != x[2] then x[4] else x[2]) > -1;
