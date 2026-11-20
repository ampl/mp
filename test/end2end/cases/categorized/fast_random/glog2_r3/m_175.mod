var x{1..6} integer >= -2 <= 0;
minimize o: (-4)*x[1] + (5)*x[2] + (3)*x[3] + (1)*x[4] + (0)*x[5] + (6)*x[6];
s.t. c1: (if x[6] == 0 then x[5] else x[5]) < -1;
s.t. c2: (x[1] != -2) ==> (x[5] < -2);
s.t. c3: (x[4] == -2) ==> (x[1] == -2) else (x[5] != -1);
s.t. c4: exists {i in 1..6} x[i] == 0;
s.t. c5: exists {i in 1..6} x[i] == 0;
