var x{1..5} integer >= 0 <= 3;
minimize o: (0)*x[1] + (-5)*x[2] + (6)*x[3] + (-1)*x[4] + (4)*x[5];
s.t. c1: exists {i in 1..5} x[i] == 2;
