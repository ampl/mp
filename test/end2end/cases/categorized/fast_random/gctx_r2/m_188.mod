var x{1..3} integer >= 0 <= 2;
minimize o: (-3)*x[1] + (0)*x[2] + (0)*x[3] + (3)*(if min(x[1], x[1]) >= 1 then 1 else 0);
s.t. c1: (x[3] != 1) ==> (x[3] == 1) else (x[3] != 1);
s.t. c2: (if min(x[1], x[1]) >= 1 then x[3] else 0) + (if x[3] == 1 then x[3] else 0) > 1;
s.t. c3: (if x[3] != 1 then x[3] else x[1]) >= 0;
s.t. c4: (min(x[1], x[1]) >= 1) and ((x[3] == 1) or not (x[3] == 1));
