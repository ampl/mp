var x{1..5} integer >= -2 <= 2;
maximize o: (0)*x[1] + (-1)*x[2] + (-4)*x[3] + (-4)*x[4] + (-5)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
