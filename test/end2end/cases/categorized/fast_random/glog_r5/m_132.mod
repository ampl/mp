var x{1..5} integer >= -2 <= 0;
minimize o: (5)*x[1] + (3)*x[2] + (5)*x[3] + (5)*x[4] + (1)*x[5];
s.t. c1: ((if x[2] >= -2 then 1 else 0) + (if x[1] >= -2 then 1 else 0)) >= 1;
s.t. c2: exists {i in 1..5} x[i] == -1;
