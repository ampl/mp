var x{1..5} integer >= 0 <= 4;
maximize o: (0)*x[1] + (4)*x[2] + (6)*x[3] + (-3)*x[4] + (2)*x[5];
s.t. c1: exists {i in 1..5} x[i] == 0;
