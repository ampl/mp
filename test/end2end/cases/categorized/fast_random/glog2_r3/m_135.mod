var x{1..4} integer >= 0 <= 3;
minimize o: (4)*x[1] + (6)*x[2] + (-5)*x[3] + (5)*x[4];
s.t. c1: (x[1] < 1) ==> (x[2] < 2) else (x[3] <= 1);
s.t. c2: (x[1] >= 0) or (x[3] == 1) or (x[1] == 0);
s.t. c3: exists {i in 1..4} x[i] == 0;
s.t. c4: max(x[1], x[4]) - min(x[3], x[4]) > 1;
s.t. c5: (x[3] == 1) ==> (x[2] != 2) else (x[2] >= 3);
