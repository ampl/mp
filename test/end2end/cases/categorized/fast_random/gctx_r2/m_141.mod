var x{1..3} integer >= 0 <= 2;
maximize o: (-2)*x[1] + (-1)*x[2] + (6)*x[3] + (-5)*(if x[1] != x[2] then 1 else 0);
s.t. c1: (max(x[1], x[1]) <= 2) or (x[1] != x[2]);
s.t. c2: (x[1] != x[2]) ==> (x[1] != x[2]) else (max(x[1], x[1]) <= 2);
s.t. c3: (x[1] != x[2]) and ((x[1] != x[2]) or not (x[1] != x[2]));
s.t. c4: (if x[1] != x[2] then x[2] else x[1]) >= 2;
