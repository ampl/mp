var x{1..4} integer >= 0 <= 1;
maximize o: (3)*x[1] + (5)*x[2] + (-2)*x[3] + (4)*x[4];
s.t. c1: atmost 3 {i in 1..4} (x[i] >= 0);
