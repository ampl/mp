var x{1..3} integer >= 0 <= 2;
maximize o: (-1)*x[1] + (-4)*x[2] + (5)*x[3];
s.t. c1: (x[2] + x[3] != 0) ==> (x[2] + x[3] != 0);
s.t. c2: (if max(x[1], x[2]) >= 2 then 1 else 0) + (if min(x[1], x[2]) == 0 then 1 else 0) + (if x[2] + x[3] != 0 then 1 else 0) >= 3;
s.t. c3: (min(x[1], x[2]) == 0) ==> (min(x[1], x[2]) == 0);
