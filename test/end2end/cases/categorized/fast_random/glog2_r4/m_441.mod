var x{1..5} integer >= 0 <= 4;
minimize o: (2)*x[1] + (5)*x[2] + (0)*x[3] + (5)*x[4] + (0)*x[5];
s.t. c1: max(x[1], x[3]) - min(x[1], x[2]) != 4;
