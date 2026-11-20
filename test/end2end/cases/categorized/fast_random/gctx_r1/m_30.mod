var x{1..3} integer >= 0 <= 3;
maximize o: (-1)*x[1] + (-4)*x[2] + (-2)*x[3];
s.t. c1: (if x[3] <= 2 then 1 else 0) + (if x[3] <= 0 then 1 else 0) + (if x[1] == 2 then 1 else 0) >= 3;
s.t. c2: (x[3] <= 0) ==> ((min(x[3], x[1]) < 1) or (x[1] == 2));
s.t. c3: (if x[3] <= 2 then x[2] else x[3]) >= 3;
s.t. c4: not ((x[3] <= 0) and (x[3] <= 0));
