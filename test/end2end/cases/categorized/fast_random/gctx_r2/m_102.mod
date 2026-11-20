var x{1..3} integer >= 0 <= 3;
maximize o: (-3)*x[1] + (6)*x[2] + (1)*x[3];
s.t. c1: (if x[1] + x[3] < 2 then x[2] else x[1]) == 0;
s.t. c2: (x[1] + x[3] < 2) and ((max(x[2], x[3]) >= 2) or not (x[1] + x[3] < 2));
s.t. c3: not (x[1] + x[3] < 2);
s.t. c4: (max(x[2], x[3]) >= 2) <==> (x[1] + x[3] < 2);
