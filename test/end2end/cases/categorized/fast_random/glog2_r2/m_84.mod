var x{1..3} integer >= -3 <= 1;
minimize o: (6)*x[1] + (4)*x[2] + (-4)*x[3];
s.t. c1: (x[3] >= -2) or (x[1] != -2) or (x[1] >= 0);
