var x{1..5} integer >= -2 <= 3;
maximize o: (-4)*x[1] + (6)*x[2] + (-1)*x[3] + (0)*x[4] + (-3)*x[5];
s.t. c1: (x[5] >= -2) ==> (x[4] <= 3);
s.t. c2: atleast 3 {i in 1..5} (x[i] >= -2);
s.t. c3: (x[5] <= 2) ==> (x[1] <= 3);
s.t. c4: (if x[5] <= -2 then x[2] else x[1]) <= 2;
