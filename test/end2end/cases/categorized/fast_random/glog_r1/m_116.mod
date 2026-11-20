var x{1..3} integer >= 0 <= 3;
minimize o: (2)*x[1] + (1)*x[2] + (-2)*x[3];
s.t. c1: abs(x[2] - x[1]) != 3;
