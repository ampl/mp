var x{1..4} integer >= 0 <= 4;
maximize o: (0)*x[1] + (6)*x[2] + (-5)*x[3] + (1)*x[4];
s.t. c1: exists {i in 1..4} x[i] == 0;
s.t. c2: ((if x[4] > 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c3: (x[1] < 3) <==> (x[2] + x[1] != 0);
s.t. c4: alldiff(x[2], x[4]);
