var x{1..5} integer >= -2 <= 1;
maximize o: (2)*x[1] + (1)*x[2] + (-5)*x[3] + (5)*x[4] + (3)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
