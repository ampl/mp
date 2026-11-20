var x{1..6} integer >= 0 <= 2;
minimize o: (2)*x[1] + (-3)*x[2] + (2)*x[3] + (3)*x[4] + (-4)*x[5] + (6)*x[6];
s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: not ((x[6] >= 1) and (x[1] > 2));
s.t. c3: (if x[3] == 2 then x[4] else x[4]) != 1;
