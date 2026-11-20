var x{1..4} integer >= 0 <= 3;
maximize o: (5)*x[1] + (5)*x[2] + (3)*x[3] + (1)*x[4];
s.t. c1: (min(x[2], x[2]) != 0) ==> (max(x[2], x[1]) > 1);
s.t. c2: (max(x[2], x[1]) > 1) ==> (max(x[2], x[1]) > 1);
s.t. c3: (if max(x[2], x[1]) > 1 then x[3] else 0) + (if min(x[2], x[2]) != 0 then x[3] else 0) >= 3;
s.t. c4: (min(x[2], x[2]) != 0) ==> (max(x[2], x[1]) > 1) else (max(x[2], x[1]) > 1);
