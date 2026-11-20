var x{1..4} integer >= 0 <= 4;
minimize o: (6)*x[1] + (3)*x[2] + (-2)*x[3] + (3)*x[4];
s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: (if x[4] != 3 then x[2] else x[1]) <= 0;
s.t. c3: (x[3] >= 2) ==> (x[1] < 0);
