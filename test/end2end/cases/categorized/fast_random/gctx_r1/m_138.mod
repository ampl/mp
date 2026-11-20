var x{1..3} integer >= -1 <= 1;
maximize o: (4)*x[1] + (3)*x[2] + (6)*x[3];
s.t. c1: (max(x[1], x[2]) < 0) ==> (max(x[1], x[2]) < 0) else (max(x[1], x[2]) < 0);
s.t. c2: (x[2] == 0) ==> ((x[2] == 0) or (max(x[1], x[2]) < 0));
s.t. c3: (if max(x[1], x[2]) < 0 then x[1] else 0) + (if x[2] == 0 then x[3] else 0) > 1;
