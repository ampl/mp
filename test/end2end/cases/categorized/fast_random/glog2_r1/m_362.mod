var x{1..5} integer >= 0 <= 3;
maximize o: (-2)*x[1] + (-2)*x[2] + (0)*x[3] + (-4)*x[4] + (2)*x[5];
s.t. c1: (x[2] == 3) and (x[2] > 2) ==> (x[2] >= 3);
s.t. c2: alldiff(x[3], x[4], x[1]);
s.t. c3: forall {i in 1..4} x[i] != x[i+1];
