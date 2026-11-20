var x{1..4} integer >= 0 <= 3;
maximize o: (1)*x[1] + (2)*x[2] + (-3)*x[3] + (0)*x[4];
s.t. c1: (x[3] <= 2) and (x[4] != 3) ==> (x[1] >= 2);
s.t. c2: (if x[1] >= 3 then x[2] else x[3]) <= 2;
