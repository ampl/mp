var x{1..4} integer >= 0 <= 3;
minimize o: (2)*x[1] + (5)*x[2] + (-3)*x[3] + (-5)*x[4];
s.t. c1: (x[4] + x[3] >= 3) ==> (min(x[3], x[2]) == 3);
s.t. c2: (if min(x[3], x[2]) == 3 then x[3] else 0) + (if x[4] + x[3] >= 3 then x[1] else 0) >= 0;
s.t. c3: not ((x[2] >= 3) and (x[2] >= 3));
