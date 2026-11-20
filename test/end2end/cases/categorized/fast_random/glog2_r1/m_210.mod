var x{1..3} integer >= 0 <= 4;
minimize o: (-5)*x[1] + (-5)*x[2] + (-3)*x[3];
s.t. c1: ((if x[3] <= 4 then 1 else 0) + (if x[3] >= 4 then 1 else 0)) >= 1;
