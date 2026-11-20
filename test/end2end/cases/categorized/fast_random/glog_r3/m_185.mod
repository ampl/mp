var x{1..4} integer >= -2 <= 2;
maximize o: (-4)*x[1] + (2)*x[2] + (5)*x[3] + (5)*x[4];
s.t. c1: (x[2] == 0) ==> (x[3] == -1) else (x[4] >= -1);
s.t. c2: (if x[1] <= 2 then x[2] else x[4]) <= 1;
