var x{1..5} integer >= 0 <= 3;
minimize o: (-5)*x[1] + (0)*x[2] + (5)*x[3] + (-1)*x[4] + (-5)*x[5];
s.t. c1: exists {i in 1..5} x[i] == 1;
s.t. c2: (x[1] >= 1) or (x[1] >= 2) or (x[1] > 2);
