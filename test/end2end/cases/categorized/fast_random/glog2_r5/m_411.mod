var x{1..6} integer >= 0 <= 4;
maximize o: (3)*x[1] + (3)*x[2] + (-4)*x[3] + (-1)*x[4] + (-3)*x[5] + (-1)*x[6];
s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: (x[6] < 0) ==> (x[3] != 0) else (x[6] <= 5);
s.t. c3: alldiff(x[5], x[1]);
