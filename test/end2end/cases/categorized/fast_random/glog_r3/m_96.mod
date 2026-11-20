var x{1..5} integer >= 0 <= 3;
minimize o: (0)*x[1] + (2)*x[2] + (-4)*x[3] + (5)*x[4] + (5)*x[5];
s.t. c1: max(x[3], x[1]) - min(x[5], x[3]) == 3;
