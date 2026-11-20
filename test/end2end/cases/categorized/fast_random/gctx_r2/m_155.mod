var x{1..3} integer >= 0 <= 2;
maximize o: (-1)*x[1] + (-4)*x[2] + (2)*x[3] + (-4)*(if x[2] == x[2] then 1 else 0);
s.t. c1: (x[2] + x[3] <= 2) <==> (x[2] + x[3] <= 2);
s.t. c2: (if min(x[2], x[1]) <= 2 then x[1] else x[3]) > 2;
s.t. c3: (if x[2] + x[3] <= 2 then x[1] else x[2]) < 0;
s.t. c4: (x[2] + x[3] <= 2) ==> ((x[2] == x[2]) or (x[2] + x[3] <= 2));
