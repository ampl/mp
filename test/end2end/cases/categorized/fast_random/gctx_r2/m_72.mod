var x{1..4} integer >= 0 <= 2;
minimize o: (-5)*x[1] + (-5)*x[2] + (-2)*x[3] + (-5)*x[4];
s.t. c1: (if min(x[3], x[4]) < 1 then x[2] else x[2]) > 2;
s.t. c2: (min(x[1], x[4]) == 0) ==> (min(x[1], x[4]) == 0) else (min(x[1], x[4]) == 0);
s.t. c3: (if min(x[1], x[4]) == 0 then x[4] else x[2]) < 2;
