var x{1..4} integer >= 0 <= 4;
minimize o: (-2)*x[1] + (1)*x[2] + (6)*x[3] + (0)*x[4];
s.t. c1: alldiff(x[3], x[4]);
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: (x[4] >= 2) or (x[4] < 1) or (x[4] > 1);
s.t. c4: (x[1] < 0) ==> (x[1] != 4) else (x[2] > 2);
