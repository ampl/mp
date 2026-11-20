var x{1..6} integer >= 0 <= 3;
minimize o: (-3)*x[1] + (5)*x[2] + (6)*x[3] + (0)*x[4] + (6)*x[5] + (5)*x[6];
s.t. c1: exists {i in 1..6} x[i] == 3;
s.t. c2: atmost 3 {i in 1..6} (x[i] >= 0);
s.t. c3: exists {i in 1..6} x[i] == 3;
s.t. c4: atmost 3 {i in 1..6} (x[i] >= 0);
