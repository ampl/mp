var x{1..5} integer >= 0 <= 1;
minimize o: (0)*x[1] + (-2)*x[2] + (4)*x[3] + (-5)*x[4] + (2)*x[5];
s.t. c1: max(x[1], x[2]) - min(x[3], x[1]) < 1;
s.t. c2: max(x[4], x[5]) - min(x[2], x[4]) > 0;
s.t. c3: (x[4] == 1) or (x[4] != 0) or (x[4] >= 0);
