var x{1..5} integer >= 0 <= 2;
maximize o: (-2)*x[1] + (6)*x[2] + (2)*x[3] + (-5)*x[4] + (2)*x[5];
s.t. c1: max(x[1], x[3]) - min(x[2], x[4]) >= 1;
s.t. c2: not ((x[2] >= 2) and (x[2] >= 0));
s.t. c3: (x[2] > 2) and (x[5] <= 1) ==> (x[5] >= 2);
s.t. c4: forall {i in 1..4} x[i] != x[i+1];
