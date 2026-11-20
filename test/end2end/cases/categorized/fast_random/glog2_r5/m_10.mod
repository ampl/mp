var x{1..6} integer >= -2 <= 2;
minimize o: (-4)*x[1] + (6)*x[2] + (-3)*x[3] + (-1)*x[4] + (-1)*x[5] + (-3)*x[6];
s.t. c1: exists {i in 1..6} x[i] == 1;
s.t. c2: ((if x[5] == 1 then 1 else 0) + (if x[6] >= 1 then 1 else 0)) >= 1;
s.t. c3: exists {i in 1..6} x[i] == -2;
s.t. c4: (x[2] <= 1) ==> (x[2] > -2);
