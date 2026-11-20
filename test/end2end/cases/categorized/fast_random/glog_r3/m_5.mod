var x{1..5} integer >= 0 <= 3;
minimize o: (0)*x[1] + (-2)*x[2] + (-2)*x[3] + (-3)*x[4] + (4)*x[5];
s.t. c1: abs(x[3] - x[1]) != 0;
