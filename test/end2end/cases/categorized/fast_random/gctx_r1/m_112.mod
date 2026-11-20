var x{1..4} integer >= 0 <= 3;
maximize o: (-2)*x[1] + (-2)*x[2] + (-2)*x[3] + (4)*x[4];
s.t. c1: (min(x[3], x[2]) != 0) or (x[3] >= x[3]);
s.t. c2: (min(x[3], x[2]) != 0) <==> (min(x[3], x[2]) != 0);
s.t. c3: (x[3] >= x[3]) ==> (min(x[3], x[2]) != 0) else (abs(x[2] - x[3]) < 1);
s.t. c4: (if min(x[3], x[2]) != 0 then x[3] else 0) + (if x[3] >= x[3] then x[4] else 0) < 0;
