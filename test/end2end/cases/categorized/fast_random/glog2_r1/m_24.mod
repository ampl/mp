var x{1..3} integer >= -2 <= 3;
minimize o: (0)*x[1] + (6)*x[2] + (5)*x[3];
s.t. c1: (x[1] > 1) ==> (x[3] == -1) else (x[2] > 0);
s.t. c2: atleast 3 {i in 1..3} (x[i] >= 1);
s.t. c3: (x[1] == 1) ==> (x[3] < 2) else (x[3] == 1);
s.t. c4: (x[2] != 3) and (x[2] >= 2) ==> (x[3] >= 3);
