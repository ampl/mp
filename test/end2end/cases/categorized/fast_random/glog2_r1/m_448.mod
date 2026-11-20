var x{1..7} integer >= 0 <= 1;
minimize o: (-4)*x[1] + (6)*x[2] + (-3)*x[3] + (5)*x[4] + (4)*x[5] + (-5)*x[6] + (-4)*x[7];
s.t. c1: ((if x[1] < 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
s.t. c2: abs(x[2] - x[2]) > 1;
s.t. c3: atmost 2 {i in 1..7} (x[i] >= 0);
