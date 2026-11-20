var x{1..3} integer >= 0 <= 2;
maximize o: (5)*x[1] + (6)*x[2] + (1)*x[3];
s.t. c1: (if x[3] >= 2 then x[2] else x[2]) == 0;
s.t. c2: (x[3] == x[1]) and ((x[3] == x[1]) or not (x[3] == x[1]));
s.t. c3: (x[3] >= 2) <==> (x[3] == x[1]);
s.t. c4: (if x[3] >= 2 then x[3] else 0) + (if x[3] >= 2 then x[1] else 0) < 0;
