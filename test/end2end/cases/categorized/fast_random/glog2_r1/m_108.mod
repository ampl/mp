var x{1..6} integer >= 0 <= 1;
minimize o: (0)*x[1] + (2)*x[2] + (-3)*x[3] + (-2)*x[4] + (6)*x[5] + (-4)*x[6];
s.t. c1: exists {i in 1..6} x[i] == 1;
s.t. c2: (x[5] != 0) ==> (x[2] != 0) else (x[6] == 1);
s.t. c3: max(x[6], x[1]) - min(x[2], x[5]) < 1;
