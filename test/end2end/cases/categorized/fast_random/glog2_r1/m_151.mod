var x{1..3} integer >= -2 <= 2;
maximize o: (3)*x[1] + (-2)*x[2] + (5)*x[3];
s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: max(x[2], x[3]) - min(x[3], x[3]) <= -2;
