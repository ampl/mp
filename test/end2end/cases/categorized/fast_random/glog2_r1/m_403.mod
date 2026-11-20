var x{1..6} integer >= -2 <= 1;
maximize o: (6)*x[1] + (-5)*x[2] + (1)*x[3] + (-1)*x[4] + (3)*x[5] + (-3)*x[6];
s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: x[5] == -2;
