var x{1..6} integer >= -2 <= 1;
maximize o: (2)*x[1] + (-1)*x[2] + (6)*x[3] + (2)*x[4] + (6)*x[5] + (4)*x[6];
s.t. c1: atleast 1 {i in 1..6} (x[i] >= -1);
s.t. c2: max(x[1], x[4]) - min(x[3], x[4]) == 0;
s.t. c3: atmost 3 {i in 1..6} (x[i] >= 0);
