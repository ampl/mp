var x{1..4} integer >= -2 <= 1;
minimize o: (5)*x[1] + (-5)*x[2] + (-4)*x[3] + (-5)*x[4];
s.t. c1: abs(x[2] - x[2]) < 2;
s.t. c2: ((if x[3] >= 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
