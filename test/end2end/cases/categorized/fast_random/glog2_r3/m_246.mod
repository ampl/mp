var x{1..6} integer >= 0 <= 2;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5] + (0)*x[6];
s.t. c1: max(x[4], x[4]) - min(x[2], x[5]) >= 0;
s.t. c2: not ((x[2] >= 2) and (x[6] >= 1));
