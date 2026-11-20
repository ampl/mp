var x{1..5} integer >= 0 <= 2;
maximize o: (6)*x[1] + (0)*x[2] + (6)*x[3] + (-2)*x[4] + (6)*x[5];
s.t. c1: (x[2] == 0) and (x[5] <= 1) ==> (x[3] >= 0);
s.t. c2: (if x[4] > 2 then x[2] else x[1]) <= 2;
s.t. c3: ((if x[4] < 2 then 1 else 0) + (if x[2] >= 2 then 1 else 0)) >= 1;
s.t. c4: atmost 2 {i in 1..5} (x[i] >= 0);
