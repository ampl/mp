var x{1..5} integer >= 0 <= 2;
maximize o: (0)*x[1] + (1)*x[2] + (2)*x[3] + (-3)*x[4] + (-1)*x[5];
s.t. c1: max(x[2], x[3]) - min(x[5], x[2]) < 0;
s.t. c2: exists {i in 1..5} x[i] == 0;
s.t. c3: (x[2] != 1) ==> (x[2] != 0) else (x[1] == 2);
s.t. c4: alldiff(x[1], x[3], x[4], x[5]);
s.t. c5: (x[4] != 1) or (x[4] >= 2) or (x[3] == 2);
