var x{1..4} integer >= 0 <= 3;
maximize o: (-3)*x[1] + (-5)*x[2] + (-3)*x[3] + (0)*x[4];
s.t. c1: (x[3] <= 2) <==> (x[3] <= 2);
s.t. c2: (if x[3] <= 2 then x[4] else x[1]) >= 1;
s.t. c3: (if x[3] <= 2 then x[3] else x[1]) < 2;
s.t. c4: (x[4] == 3) ==> (max(x[1], x[3]) < 3) else (x[3] <= 2);
