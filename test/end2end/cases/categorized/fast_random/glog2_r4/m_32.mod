var x{1..4} integer >= -2 <= 1;
maximize o: (0)*x[1] + (-3)*x[2] + (3)*x[3] + (1)*x[4];
s.t. c1: (x[4] == 0) ==> (x[3] < -2);
s.t. c2: (if x[3] == -1 then x[3] else x[3]) <= -1;
