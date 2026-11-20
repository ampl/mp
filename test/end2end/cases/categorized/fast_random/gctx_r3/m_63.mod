var x{1..4} integer >= -1 <= 2;
maximize o: (2)*x[1] + (1)*x[2] + (2)*x[3] + (4)*x[4] + (-5)*(if max(x[1], x[2]) >= 0 then 1 else 0);
s.t. c1: (max(x[1], x[2]) >= 0) ==> ((min(x[3], x[4]) < 0) or (x[1] != x[1]));
s.t. c2: (if min(x[3], x[4]) < 0 then x[3] else 0) + (if max(x[1], x[2]) >= 0 then x[1] else 0) < 0;
s.t. c3: (min(x[3], x[4]) < 0) <==> (min(x[3], x[4]) < 0);
s.t. c4: not ((min(x[2], x[1]) > 0) and (x[1] != x[1]));
