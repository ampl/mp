var x{1..4} integer >= -2 <= 0;
minimize o: (6)*x[1] + (5)*x[2] + (4)*x[3] + (-4)*x[4];
s.t. c1: exists {i in 1..4} x[i] == -1;
