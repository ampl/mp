var x{1..4} integer >= 0 <= 1;
maximize o: (-1)*x[1] + (-4)*x[2] + (3)*x[3] + (-3)*x[4];
s.t. c1: (x[2] != 0) ==> (x[4] < 1);
s.t. c2: (x[2] < 0) ==> (x[3] >= 1) else (x[3] == 1);
s.t. c3: abs(x[2] - x[3]) <= 0;
s.t. c4: (if x[1] < 1 then x[3] else x[4]) < 1;
