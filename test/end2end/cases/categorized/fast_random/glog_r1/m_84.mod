var x{1..5} integer >= -2 <= 1;
minimize o: (-2)*x[1] + (2)*x[2] + (-4)*x[3] + (-4)*x[4] + (6)*x[5];
s.t. c1: exists {i in 1..5} x[i] == 0;
