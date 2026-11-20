var x{1..4} integer >= 0 <= 2;
maximize o: (5)*x[1] + (5)*x[2] + (0)*x[3] + (-4)*x[4];
s.t. c1: (if x[2] + x[4] != 2 then 1 else 0) + (if x[2] + x[4] != 2 then 1 else 0) + (if x[2] + x[4] != 2 then 1 else 0) >= 1;
s.t. c2: (max(x[1], x[3]) <= 0) ==> (x[2] + x[4] != 2);
s.t. c3: (max(x[1], x[3]) <= 0) <==> (max(x[1], x[3]) <= 0);
s.t. c4: (x[2] + x[4] != 2) ==> ((max(x[1], x[3]) <= 0) or (max(x[1], x[3]) <= 0));
