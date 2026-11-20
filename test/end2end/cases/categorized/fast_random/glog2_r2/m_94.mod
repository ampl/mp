var x{1..5} integer >= -2 <= 2;
maximize o: (6)*x[1] + (-5)*x[2] + (6)*x[3] + (-1)*x[4] + (3)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: alldiff(x[3], x[2], x[4], x[1]);
s.t. c3: (x[3] >= 0) and (x[4] <= 1) ==> (x[5] >= 0);
