var x{1..3} integer >= -1 <= 1;
maximize o: (-5)*x[1] + (-4)*x[2] + (4)*x[3];
s.t. c1: (if x[3] + x[1] > 0 then x[1] else x[2]) <= -1;
s.t. c2: (if min(x[1], x[1]) == -1 then x[3] else 0) + (if min(x[1], x[1]) == -1 then x[2] else 0) >= 0;
s.t. c3: not (min(x[1], x[1]) == -1);
