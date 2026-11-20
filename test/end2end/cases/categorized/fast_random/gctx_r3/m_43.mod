var x{1..4} integer >= 0 <= 2;
minimize o: (-4)*x[1] + (5)*x[2] + (3)*x[3] + (0)*x[4];
s.t. c1: not ((x[3] > 0) and (x[4] != 0));
s.t. c2: (if x[4] != 0 then 1 else 0) + (if abs(x[2] - x[3]) >= 1 then 1 else 0) + (if x[3] > 0 then 1 else 0) >= 2;
s.t. c3: (x[3] > 0) or (x[3] > 0);
s.t. c4: (x[4] != 0) ==> (x[1] >= 0) else (x[1] >= 0);
