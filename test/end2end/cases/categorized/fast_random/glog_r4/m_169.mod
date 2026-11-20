var x{1..5} integer >= 0 <= 4;
minimize o: (-2)*x[1] + (4)*x[2] + (-5)*x[3] + (0)*x[4] + (-3)*x[5];
s.t. c1: abs(x[4] - x[4]) != 2;
