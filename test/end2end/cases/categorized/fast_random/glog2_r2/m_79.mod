var x{1..4} integer >= 0 <= 2;
minimize o: (4)*x[1] + (4)*x[2] + (-5)*x[3] + (3)*x[4];
s.t. c1: (x[3] < 2) ==> (x[4] < 1) else (x[4] >= 1);
s.t. c2: exists {i in 1..4} x[i] == 1;
s.t. c3: (x[3] < 2) ==> (x[4] < 1) else (x[4] >= 1);
s.t. c4: exists {i in 1..4} x[i] == 1;
