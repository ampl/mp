var x{1..5} integer >= 0 <= 3;
maximize o: (-2)*x[1] + (0)*x[2] + (-1)*x[3] + (-2)*x[4] + (2)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: (x[4] > 1) and (x[1] == 2) ==> (x[4] >= 1);
