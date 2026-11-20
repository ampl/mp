var x{1..4} integer >= 0 <= 2;
maximize o: (4)*x[1] + (1)*x[2] + (6)*x[3] + (2)*x[4];
s.t. c1: (x[3] <= 2) and (x[3] <= 2) ==> (x[1] >= 2);
s.t. c2: atmost 2 {i in 1..4} (x[i] >= 0);
s.t. c3: abs(x[2] - x[4]) > 0;
