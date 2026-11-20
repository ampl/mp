var x{1..5} integer >= 0 <= 2;
minimize o: (-3)*x[1] + (0)*x[2] + (0)*x[3] + (2)*x[4] + (1)*x[5];
s.t. c1: (x[2] != 0) ==> (x[5] == 1) else (x[5] < 1);
s.t. c2: exists {i in 1..5} x[i] == 0;
s.t. c3: not ((x[5] >= 2) and (x[1] < 2));
