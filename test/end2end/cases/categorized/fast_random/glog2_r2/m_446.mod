var x{1..4} integer >= 0 <= 2;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4];
s.t. c1: (x[1] == 1) and (x[4] < 1) ==> (x[2] >= 1);
s.t. c2: max(x[4], x[3]) - min(x[2], x[3]) < 2;
s.t. c3: exists {i in 1..4} x[i] == 0;
s.t. c4: (x[2] == 0) <==> (x[4] + x[1] <= 2);
