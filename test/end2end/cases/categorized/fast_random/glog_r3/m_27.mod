var x{1..4} integer >= -2 <= 2;
minimize o: (6)*x[1] + (1)*x[2] + (-3)*x[3] + (4)*x[4];
s.t. c1: (x[4] > -2) or (x[1] == 1) or (x[3] == 0);
