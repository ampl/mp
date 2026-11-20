var x{1..5} integer >= 0 <= 2;
maximize o: (5)*x[1] + (1)*x[2] + (-3)*x[3] + (5)*x[4] + (6)*x[5];
s.t. c1: ((if x[3] > -1 then 1 else 0) + (if x[2] >= -1 then 1 else 0)) >= 1;
s.t. c2: not ((x[4] == 2) and (x[4] >= 2));
s.t. c3: atmost 1 {i in 1..5} (x[i] >= -1);
s.t. c4: (if x[2] > -1 then x[4] else x[5]) <= 2;
