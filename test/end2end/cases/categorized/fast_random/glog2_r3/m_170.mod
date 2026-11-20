var x{1..5} integer >= 0 <= 3;
minimize o: (2)*x[1] + (6)*x[2] + (-4)*x[3] + (-3)*x[4] + (-1)*x[5];
s.t. c1: ((if x[2] > 2 then 1 else 0) + (if x[5] >= 2 then 1 else 0)) >= 1;
s.t. c2: not ((x[5] != 2) and (x[3] != 1));
s.t. c3: not ((x[5] == 3) and (x[2] != 0));
