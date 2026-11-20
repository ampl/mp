var x{1..4} integer >= -2 <= 1;
maximize o: (1)*x[1] + (-3)*x[2] + (3)*x[3] + (6)*x[4];
s.t. c1: numberof -1 in (x[1], x[2], x[3], x[4]) < 1;
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: ((if x[4] > -1 then 1 else 0) + (if x[2] >= -1 then 1 else 0)) >= 1;
s.t. c4: x[4] == -1;
