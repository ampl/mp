var x{1..5} integer >= -2 <= 1;
maximize o: (6)*x[1] + (-4)*x[2] + (0)*x[3] + (-5)*x[4] + (5)*x[5];
s.t. c1: abs(x[3] - x[1]) != 3;
s.t. c2: max(x[1], x[3]) - min(x[2], x[4]) > -2;
s.t. c3: forall {i in 1..4} x[i] != x[i+1];
