var x{1..3} integer >= -2 <= 0;
maximize o: (-1)*x[1] + (-3)*x[2] + (1)*x[3];
s.t. c1: forall {i in 1..2} x[i] != x[i+1];
