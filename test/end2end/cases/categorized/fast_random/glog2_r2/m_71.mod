var x{1..4} integer >= 0 <= 2;
maximize o: (0)*x[1] + (4)*x[2] + (1)*x[3] + (-3)*x[4];
s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: (x[2] >= 2) ==> (x[2] <= 3);
s.t. c3: (x[2] < -1) ==> (x[1] < 3) else (x[1] != 3);
s.t. c4: forall {i in 1..3} x[i] != x[i+1];
