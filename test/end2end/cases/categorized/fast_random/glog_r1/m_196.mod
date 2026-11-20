var x{1..4} integer >= 0 <= 2;
minimize o: (2)*x[1] + (-5)*x[2] + (2)*x[3] + (0)*x[4];
s.t. c1: abs(x[2] - x[3]) <= 1;
