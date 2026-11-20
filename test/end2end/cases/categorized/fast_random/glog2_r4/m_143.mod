var x{1..6} integer >= 0 <= 4;
minimize o: (5)*x[1] + (-5)*x[2] + (-5)*x[3] + (4)*x[4] + (3)*x[5] + (-3)*x[6];
s.t. c1: ((if x[4] == 0 then 1 else 0) + (if x[5] >= 0 then 1 else 0)) >= 1;
s.t. c2: x[4] == 4;
