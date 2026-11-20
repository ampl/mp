var x{1..4} integer >= 0 <= 3;
minimize o: (-4)*x[1] + (-5)*x[2] + (-4)*x[3] + (4)*x[4];
s.t. c1: (x[3] + x[4] >= 1) and ((x[3] + x[4] >= 1) or not (x[3] <= 2));
s.t. c2: (x[3] + x[4] >= 1) and ((max(x[4], x[2]) == 2) or not (max(x[4], x[2]) == 2));
s.t. c3: not ((x[3] + x[4] >= 1) and (x[3] <= 2));
s.t. c4: (if x[3] <= 2 then 1 else 0) + (if max(x[4], x[2]) == 2 then 1 else 0) + (if max(x[4], x[2]) == 2 then 1 else 0) >= 3;
