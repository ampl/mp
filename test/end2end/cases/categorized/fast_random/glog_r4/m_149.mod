var x{1..6} integer >= -2 <= 2;
minimize o: (5)*x[1] + (1)*x[2] + (-1)*x[3] + (3)*x[4] + (-3)*x[5] + (4)*x[6];
s.t. c1: (x[3] != 2) or (x[4] == 0) or (x[3] > 1);
s.t. c2: ((if x[5] <= -1 then 1 else 0) + (if x[4] >= -1 then 1 else 0)) >= 1;
s.t. c3: atmost 1 {i in 1..6} (x[i] >= -2);
