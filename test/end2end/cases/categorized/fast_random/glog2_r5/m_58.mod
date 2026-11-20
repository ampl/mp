var x{1..6} integer >= 0 <= 1;
minimize o: (4)*x[1] + (0)*x[2] + (5)*x[3] + (-4)*x[4] + (4)*x[5] + (3)*x[6];
s.t. c1: not ((x[5] == 1) and (x[5] > 1));
