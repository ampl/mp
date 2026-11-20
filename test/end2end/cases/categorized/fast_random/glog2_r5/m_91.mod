var x{1..3} integer >= 0 <= 3;
maximize o: (0)*x[1] + (2)*x[2] + (6)*x[3];
s.t. c1: (x[1] < 0) and (x[3] >= 0) ==> (x[1] >= 0);
s.t. c2: alldiff(x[1], x[3], x[2]);
s.t. c3: ((if x[1] <= 3 then 1 else 0) + (if x[3] >= 3 then 1 else 0)) >= 1;
s.t. c4: exists {i in 1..3} x[i] == 0;
