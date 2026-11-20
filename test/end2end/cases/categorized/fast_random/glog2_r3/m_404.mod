var x{1..5} integer >= -2 <= 2;
maximize o: (-1)*x[1] + (2)*x[2] + (2)*x[3] + (6)*x[4] + (-3)*x[5];
s.t. c1: (x[5] > 0) and (x[3] < 1) ==> (x[4] >= 0);
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: forall {i in 1..4} x[i] != x[i+1];
