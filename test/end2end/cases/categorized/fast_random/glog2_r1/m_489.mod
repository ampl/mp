var x{1..6} integer >= -2 <= 1;
maximize o: (6)*x[1] + (3)*x[2] + (0)*x[3] + (-5)*x[4] + (-5)*x[5] + (-1)*x[6];
s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: exists {i in 1..6} x[i] == -2;
s.t. c3: max(x[1], x[4]) - min(x[4], x[4]) >= 0;
s.t. c4: forall {i in 1..5} x[i] != x[i+1];
s.t. c5: max(x[1], x[4]) - min(x[4], x[4]) >= 0;
