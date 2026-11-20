var x{1..3} integer >= 0 <= 4;
minimize o: (0)*x[1] + (3)*x[2] + (-5)*x[3];
s.t. c1: abs(x[1] - x[3]) > 2;
