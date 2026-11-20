var x{1..4} integer >= 0 <= 1;
minimize o: (-1)*x[1] + (-4)*x[2] + (2)*x[3] + (-3)*x[4];
s.t. c1: (x[1] <= 1) ==> (x[3] > 1) else (x[4] >= 1);
s.t. c2: ((if x[4] > 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
s.t. c3: not ((x[3] >= 1) and (x[2] == 1));
