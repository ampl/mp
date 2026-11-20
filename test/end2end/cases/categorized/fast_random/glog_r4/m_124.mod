var x{1..6} integer >= 0 <= 3;
maximize o: (-1)*x[1] + (-4)*x[2] + (1)*x[3] + (-2)*x[4] + (4)*x[5] + (6)*x[6];
s.t. c1: atmost 2 {i in 1..6} (x[i] >= 2);
