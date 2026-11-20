var x{1..6} integer >= -2 <= 1;
minimize o: (3)*x[1] + (0)*x[2] + (-5)*x[3] + (-4)*x[4] + (-3)*x[5] + (-1)*x[6];
s.t. c1: (x[4] > 1) ==> (x[2] < 1) else (x[6] != -2);
s.t. c2: (if x[3] == -1 then x[1] else x[6]) >= 0;
s.t. c3: exists {i in 1..6} x[i] == -1;
s.t. c4: (x[4] > 1) ==> (x[2] < 1) else (x[6] != -2);
s.t. c5: exists {i in 1..6} x[i] == -1;
