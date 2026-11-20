var x{1..6} integer >= 0 <= 3;
maximize o: (-2)*x[1] + (-2)*x[2] + (-2)*x[3] + (-4)*x[4] + (-1)*x[5] + (1)*x[6];
s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
s.t. c3: (if x[6] > 1 then x[1] else x[3]) < 2;
s.t. c4: (x[3] <= 3) ==> (x[4] >= 0) else (x[2] > 1);
