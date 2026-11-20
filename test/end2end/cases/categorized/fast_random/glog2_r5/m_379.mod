var x{1..5} integer >= -2 <= 2;
maximize o: (5)*x[1] + (4)*x[2] + (-3)*x[3] + (2)*x[4] + (-2)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
