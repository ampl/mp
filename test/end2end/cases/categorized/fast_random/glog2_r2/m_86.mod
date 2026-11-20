var x{1..6} integer >= 0 <= 3;
maximize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5] + (0)*x[6];
s.t. c1: ((if x[5] < 2 then 1 else 0) + (if x[1] >= 2 then 1 else 0)) >= 1;
s.t. c2: (x[1] > 0) ==> (x[3] <= 2) else (x[3] < 1);
s.t. c3: exists {i in 1..6} x[i] == 0;
s.t. c4: alldiff(x[2], x[4], x[1], x[3]);
