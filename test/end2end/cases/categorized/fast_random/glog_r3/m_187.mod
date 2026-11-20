var x{1..5} integer >= 0 <= 4;
maximize o: (2)*x[1] + (5)*x[2] + (-3)*x[3] + (3)*x[4] + (-5)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
