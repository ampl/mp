var x{1..6} integer >= 0 <= 2;
maximize o: (0)*x[1] + (-4)*x[2] + (0)*x[3] + (-2)*x[4] + (6)*x[5] + (3)*x[6];
s.t. c1: exists {i in 1..6} x[i] == 2;
s.t. c2: x[5] == 1;
s.t. c3: x[1] == 1;
