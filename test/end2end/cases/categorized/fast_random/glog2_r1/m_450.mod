var x{1..6} integer >= 0 <= 3;
minimize o: (6)*x[1] + (-3)*x[2] + (-1)*x[3] + (-2)*x[4] + (-5)*x[5] + (-2)*x[6];
s.t. c1: exists {i in 1..6} x[i] == 2;
