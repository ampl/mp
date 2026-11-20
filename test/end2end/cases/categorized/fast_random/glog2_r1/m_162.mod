var x{1..4} integer >= 0 <= 2;
maximize o: (5)*x[1] + (0)*x[2] + (-4)*x[3] + (-2)*x[4];
s.t. c1: max(x[4], x[2]) - min(x[4], x[4]) < 2;
s.t. c2: atmost 2 {i in 1..4} (x[i] >= 1);
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4]) < 3;
s.t. c4: not ((x[1] >= 0) and (x[1] == 2));
