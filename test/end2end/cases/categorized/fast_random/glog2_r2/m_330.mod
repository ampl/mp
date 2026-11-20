var x{1..3} integer >= 0 <= 3;
maximize o: (-1)*x[1] + (0)*x[2] + (2)*x[3];
s.t. c1: (x[2] != 0) ==> (x[3] >= 2) else (x[2] <= 3);
s.t. c2: (x[2] == 0) <==> (x[2] + x[3] == 3);
s.t. c3: (x[1] <= 3) or (x[2] < 0) or (x[3] <= 3);
s.t. c4: (if x[2] == 2 then x[1] else x[3]) > 0;
