var x{1..4} integer >= 0 <= 3;
minimize o: (5)*x[1] + (-1)*x[2] + (2)*x[3] + (-4)*x[4];
s.t. c1: (x[1] <= 2) or (x[3] < 3) or (x[4] == 0);
s.t. c2: ((if x[2] != 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c3: abs(x[2] - x[4]) > 2;
s.t. c4: x[2] == 3;
