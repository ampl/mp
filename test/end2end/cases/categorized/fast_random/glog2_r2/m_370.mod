var x{1..6} integer >= -2 <= 0;
minimize o: (1)*x[1] + (4)*x[2] + (1)*x[3] + (-3)*x[4] + (-4)*x[5] + (0)*x[6];
s.t. c1: not ((x[2] >= -2) and (x[3] > -2));
s.t. c2: alldiff(x[3], x[5]);
