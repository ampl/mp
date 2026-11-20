var x{1..5} integer >= 0 <= 3;
maximize o: (2)*x[1] + (-2)*x[2] + (-1)*x[3] + (6)*x[4] + (0)*x[5];
s.t. c1: exists {i in 1..5} x[i] == 0;
