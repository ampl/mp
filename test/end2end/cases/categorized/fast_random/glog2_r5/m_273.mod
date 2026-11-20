var x{1..4} integer >= 0 <= 2;
minimize o: (5)*x[1] + (-2)*x[2] + (-1)*x[3] + (-4)*x[4];
s.t. c1: ((if x[4] != 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
s.t. c2: x[3] == 2;
