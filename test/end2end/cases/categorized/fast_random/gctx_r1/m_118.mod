var x{1..4} integer >= 0 <= 3;
minimize o: (-1)*x[1] + (4)*x[2] + (0)*x[3] + (2)*x[4] + (5)*(if x[3] <= 0 then 1 else 0);
s.t. c1: (if x[3] <= 0 then 1 else 0) + (if x[4] < 0 then 1 else 0) + (if x[4] < 0 then 1 else 0) >= 3;
s.t. c2: (x[3] <= 0) ==> (x[3] <= 0) else (x[4] < 0);
s.t. c3: (if x[3] <= 0 then x[2] else 0) + (if x[4] < 0 then x[2] else 0) > 1;
