var x{1..5} integer >= -2 <= 1;
maximize o: (1)*x[1] + (1)*x[2] + (3)*x[3] + (6)*x[4] + (6)*x[5];
s.t. c1: (x[1] >= 0) and (x[5] < 0) ==> (x[3] >= 0);
s.t. c2: numberof 1 in (x[1], x[2], x[3], x[4], x[5]) <= 3;
s.t. c3: (x[5] == -2) or (x[5] < 1) or (x[5] >= -1);
s.t. c4: forall {i in 1..4} x[i] != x[i+1];
