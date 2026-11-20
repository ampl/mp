var x{1..5} integer >= 0 <= 1;
maximize o: (-1)*x[1] + (1)*x[2] + (-4)*x[3] + (-1)*x[4] + (5)*x[5];
s.t. c1: max(x[2], x[3]) - min(x[1], x[4]) <= 0;
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: (x[2] == 1) ==> (x[5] <= 1);
