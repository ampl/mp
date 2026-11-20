var x{1..5} integer >= -2 <= 2;
minimize o: (6)*x[1] + (1)*x[2] + (6)*x[3] + (3)*x[4] + (-1)*x[5];
s.t. c1: exists {i in 1..5} x[i] == 3;
s.t. c2: (x[1] < -3) or (x[5] <= -3) or (x[2] < 2);
