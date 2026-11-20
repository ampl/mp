var x{1..4} integer >= 0 <= 3;
maximize o: (-1)*x[1] + (6)*x[2] + (-3)*x[3] + (6)*x[4];
s.t. c1: (if x[4] >= 0 then 1 else 0) + (if x[4] >= 0 then 1 else 0) + (if x[4] >= 0 then 1 else 0) >= 1;
s.t. c2: (x[3] + x[2] < 0) ==> ((x[3] + x[2] < 0) or (x[4] >= 0));
s.t. c3: (x[3] + x[2] < 0) <==> (x[3] + x[2] < 0);
s.t. c4: (if x[4] >= 0 then x[3] else x[3]) >= 0;
