var x{1..5} integer >= 0 <= 1;
minimize o: (1)*x[1] + (4)*x[2] + (-3)*x[3] + (5)*x[4] + (3)*x[5];
s.t. c1: (x[3] < 0) and (x[4] > 0) ==> (x[2] >= 0);
s.t. c2: max(x[3], x[1]) - min(x[4], x[1]) > 0;
s.t. c3: (x[4] < 0) ==> (x[5] == 0);
s.t. c4: exists {i in 1..5} x[i] == 0;
