var x{1..6} integer >= 0 <= 2;
maximize o: (6)*x[1] + (6)*x[2] + (5)*x[3] + (-2)*x[4] + (-5)*x[5] + (4)*x[6];
s.t. c1: (x[3] <= 2) <==> (x[6] + x[6] > 2);
s.t. c2: (if x[3] <= 0 then x[1] else x[5]) >= 0;
s.t. c3: exists {i in 1..6} x[i] == 0;
s.t. c4: forall {i in 1..5} x[i] != x[i+1];
