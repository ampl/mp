var x{1..4} integer >= 0 <= 1;
minimize o: (-3)*x[1] + (-2)*x[2] + (2)*x[3] + (1)*x[4];
s.t. c1: ((if x[1] <= 0 then 1 else 0) + (if x[4] >= 0 then 1 else 0)) >= 1;
