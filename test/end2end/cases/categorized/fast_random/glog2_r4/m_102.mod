var x{1..4} integer >= -2 <= 1;
minimize o: (6)*x[1] + (-1)*x[2] + (5)*x[3] + (4)*x[4];
s.t. c1: alldiff(x[3], x[4]);
