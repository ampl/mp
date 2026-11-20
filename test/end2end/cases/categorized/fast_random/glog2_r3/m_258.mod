var x{1..4} integer >= 0 <= 1;
maximize o: (3)*x[1] + (-1)*x[2] + (3)*x[3] + (-1)*x[4];
s.t. c1: atmost 2 {i in 1..4} (x[i] >= 1);
s.t. c2: (x[3] == 0) ==> (x[2] < 0);
s.t. c3: ((if x[1] > 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
s.t. c4: (x[4] > 1) ==> (x[4] <= 0) else (x[1] != 0);
