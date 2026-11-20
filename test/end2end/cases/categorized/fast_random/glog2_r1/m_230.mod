var x{1..3} integer >= 0 <= 3;
maximize o: (0)*x[1] + (-3)*x[2] + (-4)*x[3];
s.t. c1: atmost 1 {i in 1..3} (x[i] >= 3);
