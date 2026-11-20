var x{1..6} integer >= -2 <= 1;
minimize o: (3)*x[1] + (1)*x[2] + (1)*x[3] + (-4)*x[4] + (5)*x[5] + (-1)*x[6];
s.t. c1: ((if x[1] >= -1 then 1 else 0) + (if x[3] >= -1 then 1 else 0)) >= 1;
s.t. c2: not ((x[3] != -2) and (x[3] >= 0));
