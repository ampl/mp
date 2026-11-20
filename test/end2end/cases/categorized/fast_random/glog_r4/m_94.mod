var x{1..6} integer >= 0 <= 2;
minimize o: (6)*x[1] + (5)*x[2] + (2)*x[3] + (-2)*x[4] + (-5)*x[5] + (-5)*x[6];
s.t. c1: alldiff(x[1], x[2], x[3], x[6]);
s.t. c2: forall {i in 1..5} x[i] != x[i+1];
