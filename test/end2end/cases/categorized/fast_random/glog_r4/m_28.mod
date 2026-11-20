var x{1..4} integer >= 0 <= 4;
minimize o: (-4)*x[1] + (-1)*x[2] + (-5)*x[3] + (2)*x[4];
s.t. c1: ((if x[2] > 3 then 1 else 0) + (if x[1] >= 3 then 1 else 0)) >= 1;
