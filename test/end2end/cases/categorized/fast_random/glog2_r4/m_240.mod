var x{1..6} integer >= -2 <= 2;
maximize o: (-3)*x[1] + (-1)*x[2] + (-3)*x[3] + (-2)*x[4] + (6)*x[5] + (4)*x[6];
s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: not ((x[3] < -2) and (x[3] < 0));
