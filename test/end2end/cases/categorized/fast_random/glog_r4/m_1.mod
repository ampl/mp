var x{1..5} integer >= -2 <= 2;
minimize o: (4)*x[1] + (-3)*x[2] + (0)*x[3] + (0)*x[4] + (-3)*x[5];
s.t. c1: max(x[1], x[3]) - min(x[3], x[3]) == -2;
