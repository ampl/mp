var x{1..3} integer >= 0 <= 3;
minimize o: (1)*x[1] + (6)*x[2] + (-3)*x[3] + (-5)*(if x[1] <= x[1] then 1 else 0);
s.t. c1: (max(x[2], x[1]) == 0) ==> (x[1] <= x[1]) else (x[2] == x[3]);
s.t. c2: (max(x[2], x[1]) == 0) ==> (x[2] == x[3]);
s.t. c3: (x[1] <= x[1]) or (x[1] <= x[1]);
s.t. c4: not (x[1] <= x[1]);
