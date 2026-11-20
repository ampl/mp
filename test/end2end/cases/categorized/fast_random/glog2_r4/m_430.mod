var x{1..4} integer >= 0 <= 4;
maximize o: (-5)*x[1] + (0)*x[2] + (-4)*x[3] + (6)*x[4];
s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: not ((x[3] < 1) and (x[1] > 4));
s.t. c3: abs(x[1] - x[2]) >= 0;
