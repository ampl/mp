var x{1..6} integer >= 0 <= 3;
minimize o: (-1)*x[1] + (-5)*x[2] + (-1)*x[3] + (-3)*x[4] + (-1)*x[5] + (2)*x[6];
s.t. c1: ((if x[4] != 1 then 1 else 0) + (if x[5] >= 1 then 1 else 0)) >= 1;
s.t. c2: (if x[4] != 1 then x[4] else x[1]) >= 3;
s.t. c3: not ((x[2] < 0) and (x[1] >= 3));
