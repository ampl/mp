var x{1..4} integer >= -2 <= 1;
maximize o: (6)*x[1] + (-3)*x[2] + (4)*x[3] + (3)*x[4];
s.t. c1: (x[3] == -2) <==> (x[4] + x[1] < 1);
s.t. c2: atmost 3 {i in 1..4} (x[i] >= -2);
s.t. c3: max(x[2], x[1]) - min(x[4], x[3]) > 0;
