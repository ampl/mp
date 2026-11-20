var x{1..4} integer >= 0 <= 3;
maximize o: (5)*x[1] + (-3)*x[2] + (1)*x[3] + (4)*x[4];
s.t. c1: (if x[4] > 2 then x[2] else x[1]) == 3;
s.t. c2: (x[1] < 2) ==> (x[4] < 1) else (x[3] <= 2);
s.t. c3: x[1] == 1;
