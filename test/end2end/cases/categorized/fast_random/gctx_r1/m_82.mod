var x{1..3} integer >= 0 <= 3;
maximize o: (6)*x[1] + (4)*x[2] + (1)*x[3];
s.t. c1: (x[3] + x[2] == 2) <==> (x[3] + x[2] == 2);
s.t. c2: (if abs(x[3] - x[2]) != 1 then x[3] else x[2]) == 0;
s.t. c3: (x[3] + x[2] == 2) or (x[2] + x[1] <= 0);
s.t. c4: (if x[3] + x[2] == 2 then x[1] else 0) + (if x[3] + x[2] == 2 then x[2] else 0) < 0;
