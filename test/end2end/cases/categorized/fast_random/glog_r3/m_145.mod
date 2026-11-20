var x{1..5} integer >= -2 <= 0;
maximize o: (-1)*x[1] + (6)*x[2] + (0)*x[3] + (-3)*x[4] + (2)*x[5];
s.t. c1: (x[3] <= 0) ==> (x[1] < -2) else (x[4] < 0);
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
