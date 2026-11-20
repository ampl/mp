var x{1..5} integer >= 0 <= 3;
maximize o: (-2)*x[1] + (6)*x[2] + (0)*x[3] + (5)*x[4] + (4)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: (x[5] >= 2) ==> (x[2] < 2) else (x[3] >= 0);
s.t. c3: atmost 3 {i in 1..5} (x[i] >= 1);
s.t. c4: (x[5] <= 0) or (x[1] >= 3) or (x[5] != 3);
s.t. c5: forall {i in 1..4} x[i] != x[i+1];
s.t. c6: (x[5] <= 0) or (x[1] >= 3) or (x[5] != 3);
