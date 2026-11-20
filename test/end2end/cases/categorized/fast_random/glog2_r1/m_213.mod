var x{1..6} integer >= 0 <= 2;
minimize o: (0)*x[1] + (-4)*x[2] + (6)*x[3] + (-2)*x[4] + (-4)*x[5] + (-1)*x[6];
s.t. c1: max(x[6], x[5]) - min(x[6], x[2]) <= 2;
s.t. c2: alldiff(x[1], x[6], x[3], x[5]);
s.t. c3: x[4] == 0;
