var x{1..4} integer >= 0 <= 3;
maximize o: (-3)*x[1] + (5)*x[2] + (-4)*x[3] + (-5)*x[4];
s.t. c1: ((if x[1] >= 1 then 1 else 0) + (if x[4] >= 1 then 1 else 0)) >= 1;
s.t. c2: abs(x[4] - x[1]) == 2;
