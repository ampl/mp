var x{1..5} integer >= 0 <= 3;
minimize o: (1)*x[1] + (-5)*x[2] + (-4)*x[3] + (-5)*x[4] + (5)*x[5];
s.t. c1: not ((x[2] >= 1) and (x[4] != 3));
s.t. c2: not ((x[1] != 0) and (x[1] != 3));
s.t. c3: exists {i in 1..5} x[i] == 1;
s.t. c4: exists {i in 1..5} x[i] == 0;
