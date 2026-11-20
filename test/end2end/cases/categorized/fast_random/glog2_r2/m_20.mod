var x{1..4} integer >= 0 <= 4;
minimize o: (-1)*x[1] + (-2)*x[2] + (-1)*x[3] + (2)*x[4];
s.t. c1: exists {i in 1..4} x[i] == 4;
