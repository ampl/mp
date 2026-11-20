var x{1..4} integer >= 0 <= 2;
maximize o: (6)*x[1] + (-3)*x[2] + (5)*x[3] + (0)*x[4];
s.t. c1: (x[2] > 3) ==> (x[1] >= -1) else (x[4] > 3);
s.t. c2: (if x[4] >= -1 then x[3] else x[2]) == 3;
