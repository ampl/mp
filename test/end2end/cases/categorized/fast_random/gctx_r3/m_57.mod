var x{1..4} integer >= 0 <= 2;
minimize o: (-4)*x[1] + (-2)*x[2] + (-4)*x[3] + (2)*x[4];
s.t. c1: (if x[3] > 0 then x[4] else 0) + (if x[3] > 0 then x[2] else 0) >= 0;
s.t. c2: (x[1] <= 2) ==> (x[1] <= 2) else (x[3] > 0);
