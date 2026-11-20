var x{1..5} integer >= -2 <= 0;
minimize o: (2)*x[1] + (6)*x[2] + (6)*x[3] + (4)*x[4] + (6)*x[5];
s.t. c1: (x[1] != 0) ==> (x[5] >= -2) else (x[1] > -1);
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: exists {i in 1..5} x[i] == -1;
s.t. c4: (x[3] <= -1) and (x[4] <= 0) ==> (x[4] >= -1);
s.t. c5: x[2] == -1;
