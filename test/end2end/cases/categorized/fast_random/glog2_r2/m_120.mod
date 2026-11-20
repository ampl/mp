var x{1..3} integer >= -2 <= 0;
minimize o: (-2)*x[1] + (-3)*x[2] + (4)*x[3];
s.t. c1: ((if x[2] <= 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
