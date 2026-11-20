var x{1..3} integer >= -1 <= 1;
maximize o: (4)*x[1] + (-4)*x[2] + (-1)*x[3];
s.t. c1: (min(x[3], x[3]) > 1) <==> (x[3] == x[3]);
s.t. c2: (if x[3] == x[3] then x[1] else 0) + (if x[3] == x[3] then x[1] else 0) >= 0;
s.t. c3: (if x[3] == x[3] then x[2] else 0) + (if x[3] == x[3] then x[2] else 0) > 0;
s.t. c4: (x[3] == x[3]) ==> (min(x[3], x[3]) > 1);
