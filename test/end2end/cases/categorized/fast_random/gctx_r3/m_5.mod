var x{1..4} integer >= 0 <= 3;
maximize o: (-4)*x[1] + (3)*x[2] + (-4)*x[3] + (-5)*x[4];
s.t. c1: (if x[2] >= 0 then x[2] else x[1]) > 2;
s.t. c2: (min(x[2], x[3]) >= 3) <==> (min(x[2], x[3]) >= 3);
