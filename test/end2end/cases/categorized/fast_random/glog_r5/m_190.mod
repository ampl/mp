var x{1..4} integer >= 0 <= 2;
maximize o: (-2)*x[1] + (6)*x[2] + (-5)*x[3] + (1)*x[4];
s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: not ((x[2] == 1) and (x[1] <= 1));
