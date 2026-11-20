var x{1..4} integer >= -3 <= 2;
maximize o: (5)*x[1] + (-3)*x[2] + (6)*x[3] + (-2)*x[4];
s.t. c1: (x[4] <= 1) and (x[4] != -2) ==> (x[4] >= 1);
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
