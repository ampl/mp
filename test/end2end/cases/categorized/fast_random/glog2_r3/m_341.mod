var x{1..3} integer >= 0 <= 3;
maximize o: (-5)*x[1] + (-4)*x[2] + (4)*x[3];
s.t. c1: forall {i in 1..2} x[i] != x[i+1];
