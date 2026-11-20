var x{1..5} integer >= -3 <= 2;
maximize o: (-4)*x[1] + (2)*x[2] + (-2)*x[3] + (-3)*x[4] + (-1)*x[5];
s.t. c1: exists {i in 1..5} x[i] == 1;
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
