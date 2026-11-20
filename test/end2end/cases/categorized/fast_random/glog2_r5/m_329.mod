var x{1..6} integer >= 0 <= 4;
maximize o: (-3)*x[1] + (0)*x[2] + (3)*x[3] + (-5)*x[4] + (3)*x[5] + (-4)*x[6];
s.t. c1: exists {i in 1..6} x[i] == 1;
s.t. c2: exists {i in 1..6} x[i] == 1;
