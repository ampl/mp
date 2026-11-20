var x{1..4} integer >= 0 <= 3;
minimize o: (5)*x[1] + (-4)*x[2] + (-4)*x[3] + (4)*x[4];
s.t. c1: (x[3] <= 1) or (x[3] > 2) or (x[4] == 3);
