var x{1..5} integer >= -1 <= 2;
maximize o: (-3)*x[1] + (0)*x[2] + (2)*x[3] + (6)*x[4] + (1)*x[5];
s.t. c1: (if max(x[3], x[5]) > 2 then x[5] else x[3]) == 2;
s.t. c2: (max(x[3], x[5]) > 2) or (max(x[3], x[5]) > 2);
s.t. c3: (x[4] != x[1]) ==> (x[4] != x[1]) else (max(x[3], x[5]) > 2);
