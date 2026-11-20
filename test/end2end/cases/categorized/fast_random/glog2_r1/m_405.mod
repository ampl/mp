var x{1..6} integer >= -2 <= 2;
minimize o: (-2)*x[1] + (-5)*x[2] + (1)*x[3] + (6)*x[4] + (0)*x[5] + (4)*x[6];
s.t. c1: atmost 1 {i in 1..6} (x[i] >= -2);
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
s.t. c3: (x[3] < -1) ==> (x[1] != 1) else (x[4] < 0);
s.t. c4: (x[1] == -2) and (x[2] < 2) ==> (x[6] >= -2);
s.t. c5: (if x[1] <= 2 then x[1] else x[2]) >= 2;
