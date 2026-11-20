var x{1..3} integer >= 0 <= 2;
minimize o: (6)*x[1] + (-3)*x[2] + (-2)*x[3];
s.t. c1: not ((x[3] == 2) and (x[1] <= 0));
