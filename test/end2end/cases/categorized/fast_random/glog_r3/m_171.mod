var x{1..3} integer >= -2 <= 1;
maximize o: (5)*x[1] + (5)*x[2] + (6)*x[3];
s.t. c1: atmost 3 {i in 1..3} (x[i] >= -1);
s.t. c2: ((if x[2] == 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
s.t. c3: (x[2] < -1) or (x[1] != 1) or (x[2] <= 0);
