var x{1..3} integer >= 0 <= 2;
maximize o: (-1)*x[1] + (4)*x[2] + (-3)*x[3];
s.t. c1: (x[1] == 2) <==> (x[3] <= 1);
s.t. c2: (max(x[2], x[2]) >= 1) <==> (max(x[2], x[2]) >= 1);
s.t. c3: (if max(x[2], x[2]) >= 1 then x[2] else 0) + (if x[1] == 2 then x[1] else 0) >= 2;
s.t. c4: (x[1] == 2) <==> (x[3] <= 1);
