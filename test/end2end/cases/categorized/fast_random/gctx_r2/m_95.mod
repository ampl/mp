var x{1..3} integer >= 0 <= 2;
maximize o: (2)*x[1] + (1)*x[2] + (-3)*x[3];
s.t. c1: (if min(x[3], x[1]) >= 2 then 1 else 0) + (if x[1] + x[2] < 2 then 1 else 0) + (if x[2] >= x[1] then 1 else 0) >= 2;
s.t. c2: (x[1] + x[2] < 2) ==> (x[1] + x[2] < 2);
