var x{1..5} integer >= 0 <= 3;
minimize o: (6)*x[1] + (-1)*x[2] + (-3)*x[3] + (-3)*x[4] + (-5)*x[5];
s.t. c1: abs(x[1] - x[1]) >= 0;
