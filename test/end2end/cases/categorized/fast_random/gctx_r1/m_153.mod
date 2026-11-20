var x{1..4} integer >= -1 <= 2;
minimize o: (0)*x[1] + (-2)*x[2] + (3)*x[3] + (3)*x[4] + (-5)*(if x[3] + x[3] <= 1 then 1 else 0);
s.t. c1: (x[4] != -1) ==> (x[4] != -1) else (x[3] + x[3] <= 1);
s.t. c2: (if x[4] != -1 then x[3] else 0) + (if x[3] + x[3] <= 1 then x[1] else 0) < -1;
s.t. c3: not (x[3] + x[3] <= 1);
