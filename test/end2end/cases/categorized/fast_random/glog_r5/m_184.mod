var x{1..5} integer >= 0 <= 3;
maximize o: (2)*x[1] + (-3)*x[2] + (1)*x[3] + (4)*x[4] + (3)*x[5];
s.t. c1: (x[2] >= 2) ==> (x[5] >= 1);
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: (x[2] <= 1) and (x[1] >= 3) ==> (x[5] >= 1);
s.t. c4: forall {i in 1..4} x[i] != x[i+1];
