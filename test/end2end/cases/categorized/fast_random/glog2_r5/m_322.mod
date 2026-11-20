var x{1..3} integer >= 0 <= 3;
maximize o: (-4)*x[1] + (6)*x[2] + (0)*x[3];
s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: max(x[1], x[1]) - min(x[1], x[1]) <= 3;
