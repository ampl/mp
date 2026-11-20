var x{1..6} integer >= -2 <= 1;
minimize o: (-2)*x[1] + (3)*x[2] + (6)*x[3] + (1)*x[4] + (-3)*x[5] + (-5)*x[6];
s.t. c1: exists {i in 1..6} x[i] == 1;
s.t. c2: x[3] == -2;
