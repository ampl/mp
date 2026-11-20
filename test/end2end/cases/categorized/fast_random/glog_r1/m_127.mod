var x{1..6} integer >= 0 <= 2;
minimize o: (0)*x[1] + (6)*x[2] + (-3)*x[3] + (-5)*x[4] + (2)*x[5] + (-1)*x[6];
s.t. c1: exists {i in 1..6} x[i] == 1;
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
s.t. c3: (if x[5] < 2 then x[4] else x[1]) >= 2;
s.t. c4: (x[1] != 1) and (x[3] != 1) ==> (x[4] >= 1);
