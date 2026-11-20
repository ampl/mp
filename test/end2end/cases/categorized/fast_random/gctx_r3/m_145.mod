var x{1..3} integer >= -1 <= 2;
maximize o: (-5)*x[1] + (-3)*x[2] + (-5)*x[3];
s.t. c1: (if max(x[3], x[1]) > 1 then x[2] else 0) + (if min(x[2], x[2]) >= 2 then x[1] else 0) <= 2;
s.t. c2: (max(x[3], x[1]) > 1) <==> (max(x[3], x[1]) > 1);
s.t. c3: (max(x[2], x[2]) < -1) or (max(x[3], x[1]) > 1);
s.t. c4: not (min(x[3], x[3]) >= 0);
