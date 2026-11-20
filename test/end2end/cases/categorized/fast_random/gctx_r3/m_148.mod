var x{1..4} integer >= 0 <= 3;
minimize o: (3)*x[1] + (4)*x[2] + (-5)*x[3] + (-2)*x[4];
s.t. c1: (min(x[2], x[4]) < 1) and ((x[1] != 0) or not (x[1] != 0));
s.t. c2: (min(x[3], x[2]) < 0) or (x[1] != 0);
s.t. c3: (if min(x[3], x[2]) < 0 then x[4] else 0) + (if min(x[3], x[2]) < 0 then x[1] else 0) < 2;
s.t. c4: (min(x[3], x[2]) < 0) ==> (x[1] != 0);
