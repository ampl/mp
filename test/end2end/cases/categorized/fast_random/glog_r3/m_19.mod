var x{1..4} integer >= 0 <= 2;
maximize o: (1)*x[1] + (5)*x[2] + (2)*x[3] + (4)*x[4];
s.t. c1: (x[2] >= 0) <==> (x[2] + x[1] < 2);
s.t. c2: (x[4] <= 2) or (x[1] >= 2) or (x[4] > 2);
s.t. c3: alldiff(x[2], x[4]);
s.t. c4: atmost 1 {i in 1..4} (x[i] >= 1);
