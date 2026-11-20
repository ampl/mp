var x{1..3} integer >= 0 <= 4;
maximize o: (0)*x[1] + (4)*x[2] + (-5)*x[3];
s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
