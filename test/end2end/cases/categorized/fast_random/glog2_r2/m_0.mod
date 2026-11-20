var x{1..5} integer >= 0 <= 2;
maximize o: (2)*x[1] + (5)*x[2] + (-3)*x[3] + (0)*x[4] + (-3)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: alldiff(x[2], x[3], x[1]);
s.t. c3: not ((x[3] > 2) and (x[5] < 1));
s.t. c4: (x[4] <= 0) ==> (x[1] <= 1);
