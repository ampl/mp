var x{1..6} integer >= 0 <= 3;
minimize o: (4)*x[1] + (-3)*x[2] + (0)*x[3] + (5)*x[4] + (3)*x[5] + (-4)*x[6];
s.t. c1: abs(x[4] - x[3]) < 2;
s.t. c2: alldiff(x[2], x[4], x[5], x[1]);
