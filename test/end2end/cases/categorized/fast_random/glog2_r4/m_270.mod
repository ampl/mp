var x{1..5} integer >= -2 <= 1;
minimize o: (-2)*x[1] + (-5)*x[2] + (5)*x[3] + (2)*x[4] + (-2)*x[5];
s.t. c1: abs(x[4] - x[3]) > 2;
s.t. c2: ((if x[2] != 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
s.t. c3: abs(x[2] - x[3]) <= 2;
s.t. c4: not ((x[4] >= 1) and (x[5] != -1));
