var x{1..6} integer >= -2 <= 2;
minimize o: (-1)*x[1] + (-5)*x[2] + (6)*x[3] + (-2)*x[4] + (1)*x[5] + (4)*x[6];
s.t. c1: max(x[4], x[6]) - min(x[2], x[6]) > 0;
s.t. c2: exists {i in 1..6} x[i] == 2;
s.t. c3: (x[4] > 0) ==> (x[3] < -1) else (x[5] < 0);
s.t. c4: (if x[4] != -2 then x[6] else x[6]) != 2;
s.t. c5: forall {i in 1..5} x[i] != x[i+1];
s.t. c6: (if x[4] != 1 then x[1] else x[6]) > -2;
