var x{1..5} integer >= -2 <= 3;
maximize o: (-3)*x[1] + (6)*x[2] + (1)*x[3] + (0)*x[4] + (-4)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: (if x[3] > -1 then x[5] else x[1]) < 1;
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) >= 1;
s.t. c4: (x[4] >= -1) ==> (x[4] != -1) else (x[1] > 2);
