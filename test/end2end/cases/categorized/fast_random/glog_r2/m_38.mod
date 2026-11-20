var x{1..6} integer >= 0 <= 2;
maximize o: (-5)*x[1] + (6)*x[2] + (5)*x[3] + (-2)*x[4] + (3)*x[5] + (5)*x[6];
s.t. c1: exists {i in 1..6} x[i] == 0;
s.t. c2: ((if x[2] < 1 then 1 else 0) + (if x[5] >= 1 then 1 else 0)) >= 1;
s.t. c3: exists {i in 1..6} x[i] == 1;
s.t. c4: (x[2] > 1) or (x[5] >= 1) or (x[6] < 1);
