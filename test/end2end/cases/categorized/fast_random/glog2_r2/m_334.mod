var x{1..5} integer >= -3 <= 2;
maximize o: (5)*x[1] + (-1)*x[2] + (-3)*x[3] + (2)*x[4] + (1)*x[5];
s.t. c1: (x[4] != -2) and (x[2] <= -1) ==> (x[5] >= -2);
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
