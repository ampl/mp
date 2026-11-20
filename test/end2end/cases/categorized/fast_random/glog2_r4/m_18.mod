var x{1..7} integer >= 0 <= 1;
maximize o: (1)*x[1] + (-5)*x[2] + (-4)*x[3] + (-2)*x[4] + (5)*x[5] + (-3)*x[6] + (0)*x[7];
s.t. c1: (x[1] == 1) or (x[6] > 0) or (x[1] < 1);
s.t. c2: ((if x[7] <= 1 then 1 else 0) + (if x[7] >= 1 then 1 else 0)) >= 1;
