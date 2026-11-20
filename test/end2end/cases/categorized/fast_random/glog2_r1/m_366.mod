var x{1..4} integer >= -2 <= 2;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4];
s.t. c1: exists {i in 1..4} x[i] == 0;
