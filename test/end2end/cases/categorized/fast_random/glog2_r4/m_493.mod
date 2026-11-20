var x{1..5} integer >= 0 <= 3;
minimize o: (6)*x[1] + (5)*x[2] + (-4)*x[3] + (-1)*x[4] + (6)*x[5];
s.t. c1: exists {i in 1..5} x[i] == 2;
s.t. c2: x[2] == 2;
