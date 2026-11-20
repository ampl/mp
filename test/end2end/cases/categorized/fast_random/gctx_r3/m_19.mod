var x{1..5} integer >= 0 <= 3;
minimize o: (-4)*x[1] + (-5)*x[2] + (-3)*x[3] + (2)*x[4] + (-4)*x[5] + (-2)*(if x[4] >= 0 then 1 else 0);
s.t. c1: (x[4] >= 0) and ((max(x[1], x[2]) <= 2) or not (x[2] == 3));
s.t. c2: (max(x[1], x[2]) <= 2) ==> (x[2] == 3);
