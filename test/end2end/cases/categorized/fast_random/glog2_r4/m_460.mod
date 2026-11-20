var x{1..3} integer >= -2 <= 1;
minimize o: (2)*x[1] + (5)*x[2] + (2)*x[3];
s.t. c1: abs(x[3] - x[2]) != 0;
