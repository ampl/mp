var x{1..4} integer >= 0 <= 3;
maximize o: (-4)*x[1] + (0)*x[2] + (-1)*x[3] + (3)*x[4];
s.t. c1: alldiff(x[2], x[4]);
s.t. c2: atmost 2 {i in 1..4} (x[i] >= 1);
s.t. c3: (x[3] < 2) ==> (x[3] >= 0);
s.t. c4: ((if x[2] == 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
