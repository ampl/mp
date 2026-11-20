var x{1..4} integer >= -2 <= 0;
maximize o: (-666)*x[1] + (34)*x[2] + (-999)*x[3] + (3000)*x[4];
s.t. c1: (x[3] < -2) or (x[4] == -2) or (x[1] == -1);
s.t. c2: atmost 1 {i in 1..4} (x[i] >= -2);
