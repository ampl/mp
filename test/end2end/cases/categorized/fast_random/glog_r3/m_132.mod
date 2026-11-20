var x{1..4} integer >= 0 <= 3;
maximize o: (3)*x[1] + (4)*x[2] + (-2)*x[3] + (-2)*x[4];
s.t. c1: atmost 1 {i in 1..4} (x[i] >= 2);
s.t. c2: abs(x[1] - x[3]) >= 1;
