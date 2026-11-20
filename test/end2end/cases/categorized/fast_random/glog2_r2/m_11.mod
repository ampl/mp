var x{1..5} integer >= 0 <= 4;
minimize o: (0)*x[1] + (-4)*x[2] + (2)*x[3] + (1)*x[4] + (-5)*x[5];
s.t. c1: alldiff(x[1], x[4]);
s.t. c2: (x[3] < 0) or (x[4] >= 0) or (x[1] < 4);
s.t. c3: exists {i in 1..5} x[i] == -1;
