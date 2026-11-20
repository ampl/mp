var x{1..5} integer >= -2 <= 2;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5];
s.t. c1: exists {i in 1..5} x[i] == 2;
s.t. c2: (x[4] != 0) and (x[3] >= 2) ==> (x[5] >= 0);
s.t. c3: (x[3] > 2) and (x[4] != -1) ==> (x[2] >= 2);
s.t. c4: exists {i in 1..5} x[i] == -2;
