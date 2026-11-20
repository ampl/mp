var x{1..5} integer >= 0 <= 3;
minimize o: (-1)*x[1] + (5)*x[2] + (4)*x[3] + (-3)*x[4] + (6)*x[5];
s.t. c1: exists {i in 1..5} x[i] == 3;
