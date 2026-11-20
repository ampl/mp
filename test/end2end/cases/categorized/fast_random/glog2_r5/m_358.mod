var x{1..4} integer >= 0 <= 1;
maximize o: (1)*x[1] + (-5)*x[2] + (-3)*x[3] + (4)*x[4];
s.t. c1: forall {i in 1..3} x[i] != x[i+1];
