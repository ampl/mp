var x{1..3} integer >= -1 <= 1;
maximize o: (6)*x[1] + (-2)*x[2] + (6)*x[3];
s.t. c1: (min(x[3], x[2]) <= -1) and ((x[3] == x[2]) or not (x[3] == x[2]));
s.t. c2: (if min(x[3], x[2]) <= -1 then x[3] else 0) + (if min(x[3], x[2]) <= -1 then x[1] else 0) >= 1;
s.t. c3: not (min(x[3], x[2]) <= -1);
