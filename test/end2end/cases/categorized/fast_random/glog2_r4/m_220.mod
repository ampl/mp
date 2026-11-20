var x{1..3} integer >= 0 <= 2;
minimize o: (1)*x[1] + (6)*x[2] + (-1)*x[3];
s.t. c1: exists {i in 1..3} x[i] == 1;
s.t. c2: max(x[3], x[2]) - min(x[1], x[3]) > 2;
