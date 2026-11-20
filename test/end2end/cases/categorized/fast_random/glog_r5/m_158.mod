var x{1..6} integer >= -2 <= 2;
minimize o: (1)*x[1] + (-5)*x[2] + (6)*x[3] + (-2)*x[4] + (6)*x[5] + (1)*x[6];
s.t. c1: (if x[2] == -1 then x[3] else x[5]) < -2;
s.t. c2: (x[6] == 2) ==> (x[5] > 0) else (x[6] > 2);
s.t. c3: exists {i in 1..6} x[i] == -1;
s.t. c4: max(x[3], x[3]) - min(x[4], x[2]) > -1;
