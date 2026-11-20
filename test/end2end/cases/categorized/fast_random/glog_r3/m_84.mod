var x{1..6} integer >= 0 <= 3;
minimize o: (-1)*x[1] + (0)*x[2] + (6)*x[3] + (-2)*x[4] + (5)*x[5] + (3)*x[6];
s.t. c1: (x[3] >= 2) and (x[5] == 2) ==> (x[2] >= 2);
s.t. c2: exists {i in 1..6} x[i] == 3;
