var x{1..4} integer >= 0 <= 2;
minimize o: (0)*x[1] + (-4)*x[2] + (4)*x[3] + (2)*x[4];
s.t. c1: not ((x[4] == 1) and (x[4] < 0));
