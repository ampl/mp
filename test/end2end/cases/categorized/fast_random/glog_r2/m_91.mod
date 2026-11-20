var x{1..5} integer >= -2 <= 2;
minimize o: (-1)*x[1] + (4)*x[2] + (1)*x[3] + (5)*x[4] + (3)*x[5];
s.t. c1: exists {i in 1..5} x[i] == 0;
