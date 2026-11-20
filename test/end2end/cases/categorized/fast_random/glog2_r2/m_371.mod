var x{1..3} integer >= 0 <= 4;
maximize o: (2)*x[1] + (-4)*x[2] + (6)*x[3];
s.t. c1: atmost 3 {i in 1..3} (x[i] >= 4);
s.t. c2: abs(x[1] - x[1]) < 1;
