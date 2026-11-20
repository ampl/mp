var x{1..3} integer >= 0 <= 2;
minimize o: (2)*x[1] + (-2)*x[2] + (0)*x[3];
s.t. c1: (x[2] >= 2) or (x[2] >= 1) or (x[1] <= 1);
