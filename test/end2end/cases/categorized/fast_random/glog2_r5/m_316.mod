var x{1..3} integer >= 0 <= 2;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3];
s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: (x[2] > 0) and (x[2] != 0) ==> (x[2] >= 0);
s.t. c3: (x[3] == 0) ==> (x[2] > 2) else (x[1] <= 1);
s.t. c4: abs(x[1] - x[1]) == 2;
