var x{1..6} integer >= 0 <= 2;
minimize o: (3)*x[1] + (-5)*x[2] + (6)*x[3] + (-5)*x[4] + (-2)*x[5] + (6)*x[6];
s.t. c1: alldiff(x[1], x[2], x[3]);
s.t. c2: (x[3] > 2) and (x[4] <= 2) ==> (x[3] >= 2);
s.t. c3: exists {i in 1..6} x[i] == 1;
s.t. c4: (x[2] <= 0) and (x[5] == 0) ==> (x[2] >= 0);
